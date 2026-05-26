// CurlyArrows Base JavaScript - Flask Version (UPDATED to auto-load exercise JSON into sketcher)
var debug = false;

/* global ChemDoodle, sketcher */

// ====================================================================
//  Exercise dropdown / URL param -> load JSON into ChemDoodle Sketcher
//  - CSV (semicolon) mapping: /static/exe_img/json.csv
//  - Robust JSON normalizer (Excel doubled quotes, base64, URL/HTML-escaped)
//  - Direct API load (contentFrom/molFrom) with fallback to admin UI importer
//  - Optional left-align molecules + shapes to avoid off-canvas/center placement
//  - Wait for sketcher readiness BEFORE loading
// ====================================================================

// ---------- CONFIG / DIAGNOSTICS ----------
function logStep(step, detail) {
  if (!debug) return;
  try { console.log(`[CHEM] ${step}`, detail ?? ''); } catch {}
}
window.addEventListener('error', (e) => {
  console.error('[CHEM] window error', e?.error || e);
});

// ---------- Canvas visibility helpers (optional; prevents center→left flicker) ----------
function getCanvasEl() { return document.getElementById('sketcher'); }
function hideCanvas()  { const el = getCanvasEl(); if (el) el.classList.add('is-loading'); }
function showCanvas()  { const el = getCanvasEl(); if (el) el.classList.remove('is-loading'); }

// Initialize variables
var historyData = [];
var dateTime = new Date();

// Get URL parameters
const urlParams = new URLSearchParams(window.location.search);
const exerciseParam = urlParams.get('exercise'); // used when you navigate back from another tab

// ====================================================================
//          CSV (semicolon) -> JSON -> ChemDoodle loader
// ====================================================================

// --- parsing helpers (normalize key & value) ---
function stripBOM(s){ return s && s.charCodeAt(0) === 0xFEFF ? s.slice(1) : s; }
function stripOuterQuotes(s){
  if (!s) return s;
  s = s.trim();
  const q = s[0];
  if ((q === '"' || q === "'") && s[s.length-1] === q) return s.slice(1, -1);
  return s;
}
function normalizeKey(raw){
  let k = stripBOM(raw || ''); k = k.trim(); k = stripOuterQuotes(k); return k.trim();
}
function normalizeVal(raw){
  let v = stripBOM(raw || ''); v = stripOuterQuotes(v); return v.trim();
}

let EX_JSON_MAP = new Map();

function parseCsvToMap(text) {
  const map = new Map();
  text = stripBOM(text || '');
  const lines = text.split(/\r?\n/).filter(l => l.trim().length);
  for (const line of lines) {
    const idx = line.indexOf(';');
    if (idx === -1) continue;
    const rawKey = line.slice(0, idx);
    const rawVal = line.slice(idx + 1);
    const key = normalizeKey(rawKey);
    const val = normalizeVal(rawVal);
    if (!key) continue;
    map.set(key, val);
  }
  return map;
}

async function preloadExerciseJsonMap() {
  const url = '/static/exe_img/json.csv';
  logStep('fetch csv', url);
  const res = await fetch(url, { cache: 'no-store' });
  if (!res.ok) throw new Error(`CSV fetch failed: ${res.status} ${res.statusText}`);
  const text = await res.text();
  EX_JSON_MAP = parseCsvToMap(text);
  logStep('csv parsed keys', [...EX_JSON_MAP.keys()]);
}

// --- robust value normalization (accept many encodings + Excel fix) ---
function unescapeHtmlEntities(s) {
  return s
    .replace(/&quot;/g, '"').replace(/&#34;/g, '"').replace(/&#x22;/gi, '"')
    .replace(/&amp;/g, '&').replace(/&#38;/g, '&')
    .replace(/&lt;/g, '<').replace(/&#60;/g, '<')
    .replace(/&gt;/g, '>').replace(/&#62;/g, '>')
    .replace(/&#123;/g, '{').replace(/&#125;/g, '}')
    .replace(/&#91;/g, '[').replace(/&#93;/g, ']');
}
function tryParseOnce(s) { try { JSON.parse(s); return true; } catch { return false; } }
function base64urlToStd(b64url) {
  let s = b64url.replace(/-/g, '+').replace(/_/g, '/');
  const pad = s.length % 4; if (pad) s += '='.repeat(4 - pad);
  return s;
}
function tryAtobVariants(s) {
  try { return { ok: true, value: atob(s) }; } catch {}
  try { return { ok: true, value: atob(base64urlToStd(s)) }; } catch {}
  return { ok: false };
}
function maybeDecodeURIComponent(s) {
  if (/%[0-9A-Fa-f]{2}/.test(s) || /%7B|%7D|%22/i.test(s)) {
    try { return { ok: true, value: decodeURIComponent(s) }; } catch {}
  }
  return { ok: false };
}
function collapseExcelDoubledQuotes(s) { return s.includes('""') ? s.replace(/""/g, '"') : s; }

function normalizeJsonString(rawInput) {
  let s = stripBOM(String(rawInput || '')).trim();

  const candidates = [];

  // raw & excel-dq on raw
  candidates.push({ label: 'raw', s });
  const rawExcel = collapseExcelDoubledQuotes(s);
  if (rawExcel !== s) candidates.push({ label: 'raw(excel-dq)', s: rawExcel });

  // strip outer quotes, unescape \" and doubled quotes
  let deq = stripOuterQuotes(s);
  if (deq !== s) {
    deq = deq.replace(/\\"/g, '"').replace(/""/g, '"').replace(/\\n/g, '\n').replace(/\\t/g, '\t');
    candidates.push({ label: 'dequoted', s: deq });
    const deqExcel = collapseExcelDoubledQuotes(deq);
    if (deqExcel !== deq) candidates.push({ label: 'dequoted(excel-dq)', s: deqExcel });
  }

  // html-unescaped variants
  const html1 = unescapeHtmlEntities(s);
  if (html1 !== s) {
    candidates.push({ label: 'html-unescaped', s: html1 });
    const html1Excel = collapseExcelDoubledQuotes(html1);
    if (html1Excel !== html1) candidates.push({ label: 'html-unescaped(excel-dq)', s: html1Excel });
  }
  if (deq && deq !== s) {
    const html2 = unescapeHtmlEntities(deq);
    if (html2 !== deq) {
      candidates.push({ label: 'dequoted+html', s: html2 });
      const html2Excel = collapseExcelDoubledQuotes(html2);
      if (html2Excel !== html2) candidates.push({ label: 'dequoted+html(excel-dq)', s: html2Excel });
    }
  }

  // URL-decoding
  const url1 = maybeDecodeURIComponent(s);
  if (url1.ok) {
    candidates.push({ label: 'url-decoded', s: url1.value });
    const url1Excel = collapseExcelDoubledQuotes(url1.value);
    if (url1Excel !== url1.value) candidates.push({ label: 'url-decoded(excel-dq)', s: url1Excel });
  }
  if (deq && deq !== s) {
    const url2 = maybeDecodeURIComponent(deq);
    if (url2.ok) {
      candidates.push({ label: 'dequoted+url', s: url2.value });
      const url2Excel = collapseExcelDoubledQuotes(url2.value);
      if (url2Excel !== url2.value) candidates.push({ label: 'dequoted+url(excel-dq)', s: url2Excel });
    }
  }

  // base64 -> text (and transforms on decoded text)
  const b64a = tryAtobVariants(s);
  if (b64a.ok) {
    candidates.push({ label: 'base64->text', s: b64a.value });

    const b64Excel = collapseExcelDoubledQuotes(b64a.value);
    if (b64Excel !== b64a.value) candidates.push({ label: 'base64->text(excel-dq)', s: b64Excel });

    const b64d = stripOuterQuotes(b64a.value);
    if (b64d !== b64a.value) {
      candidates.push({ label: 'base64->text(dequoted)', s: b64d });
      const b64dExcel = collapseExcelDoubledQuotes(b64d);
      if (b64dExcel !== b64d) candidates.push({ label: 'base64->text(dequoted,excel-dq)', s: b64dExcel });
    }

    const b64h = unescapeHtmlEntities(b64a.value);
    if (b64h !== b64a.value) {
      candidates.push({ label: 'base64->text(html)', s: b64h });
      const b64hExcel = collapseExcelDoubledQuotes(b64h);
      if (b64hExcel !== b64h) candidates.push({ label: 'base64->text(html,excel-dq)', s: b64hExcel });
    }

    const b64u = maybeDecodeURIComponent(b64a.value);
    if (b64u.ok) {
      candidates.push({ label: 'base64->text(url)', s: b64u.value });
      const b64uExcel = collapseExcelDoubledQuotes(b64u.value);
      if (b64uExcel !== b64u.value) candidates.push({ label: 'base64->text(url,excel-dq)', s: b64uExcel });
    }
  }

  // Attempt parses – only return strings that JSON.parse() accepts
  for (const cand of candidates) {
    if (tryParseOnce(cand.s)) {
      // double-encoded JSON string? parse twice
      let parsed;
      try { parsed = JSON.parse(cand.s); } catch {}
      if (typeof parsed === 'string' && tryParseOnce(parsed)) {
        logStep('normalizeJsonString success', cand.label + ' (double-parse)');
        return parsed;
      }
      logStep('normalizeJsonString success', cand.label);
      return cand.s;
    }
  }

  // Last attempt: force collapse of doubled quotes
  if (s.includes('""')) {
    const forced = s.replace(/""/g, '"');
    if (tryParseOnce(forced)) {
      logStep('normalizeJsonString success', 'forced excel-dq');
      return forced;
    }
  }

  throw new Error('CSV value is neither valid JSON nor base64-encoded JSON.');
}

// ====================================================================
//            OPTIONAL LEFT-ALIGN HELPERS (shift after loading)
// ====================================================================

function findMinXInObject(obj, currentMin = Infinity, seen = new Set()) {
  if (!obj || typeof obj !== 'object' || seen.has(obj)) return currentMin;
  seen.add(obj);

  if (Array.isArray(obj)) {
    for (const item of obj) currentMin = findMinXInObject(item, currentMin, seen);
    return currentMin;
  }

  for (const key in obj) {
    if (!Object.prototype.hasOwnProperty.call(obj, key)) continue;
    const val = obj[key];
    if (key === 'x' && typeof val === 'number' && isFinite(val)) {
      if (val < currentMin) currentMin = val;
    } else if (val && typeof val === 'object') {
      currentMin = findMinXInObject(val, currentMin, seen);
    }
  }
  return currentMin;
}

function shiftAllXInObject(obj, dx, seen = new Set()) {
  if (!obj || typeof obj !== 'object' || seen.has(obj)) return;
  seen.add(obj);

  if (Array.isArray(obj)) {
    for (const item of obj) shiftAllXInObject(item, dx, seen);
    return;
  }

  for (const key in obj) {
    if (!Object.prototype.hasOwnProperty.call(obj, key)) continue;
    const val = obj[key];
    if (key === 'x' && typeof val === 'number' && isFinite(val)) {
      obj[key] = val + dx;
    } else if (val && typeof val === 'object') {
      shiftAllXInObject(val, dx, seen);
    }
  }
}

function alignMoleculesLeft(sketcherInst, paddingPx = 50) {
  try {
    const mols = sketcherInst.molecules || [];
    const shapes = sketcherInst.shapes || [];
    if (!mols.length && !shapes.length) return;

    let minX = Infinity;
    for (const mol of mols) {
      if (!mol?.atoms) continue;
      for (const atom of mol.atoms) {
        if (typeof atom.x === 'number' && atom.x < minX) minX = atom.x;
      }
    }
    for (const sh of shapes) minX = findMinXInObject(sh, minX);

    if (!isFinite(minX)) return;
    const dx = paddingPx - minX;

    for (const mol of mols) {
      if (mol?.atoms) for (const atom of mol.atoms) if (typeof atom.x === 'number') atom.x += dx;
      if (Array.isArray(mol?.bonds)) {
        for (const b of mol.bonds) {
          if (typeof b.cp1x === 'number') b.cp1x += dx;
          if (typeof b.cp2x === 'number') b.cp2x += dx;
        }
      }
    }
    for (const sh of shapes) shiftAllXInObject(sh, dx);

    sketcherInst.repaint();
    logStep('aligned left', { padding: paddingPx, appliedShift: dx.toFixed(2) });
  } catch (e) {
    console.error('[CHEM] alignMoleculesLeft failed', e);
  }
}

// ====================================================================
//     ADMIN UI IMPORTER FALLBACK (if direct API fails)
// ====================================================================

function waitForEl(sel, timeoutMs = 8000) {
  return new Promise((resolve, reject) => {
    const el = document.querySelector(sel);
    if (el) return resolve(el);
    const obs = new MutationObserver(() => {
      const el2 = document.querySelector(sel);
      if (el2) { obs.disconnect(); resolve(el2); }
    });
    obs.observe(document.documentElement, { childList: true, subtree: true });
    setTimeout(() => { obs.disconnect(); reject(new Error(`Timeout waiting for ${sel}`)); }, timeoutMs);
  });
}

async function importViaAdminUI(jsonString) {
  const textEl = await waitForEl('#sketcher_open_text');
  const btnEl  = await waitForEl('#sketcher_open_load');
  textEl.value = jsonString;
  logStep('admin textarea set', jsonString.slice(0, 200) + (jsonString.length > 200 ? '…' : ''));
  btnEl.click();
  logStep('admin Open clicked');

  setTimeout(() => {
    try { alignMoleculesLeft(window.sketcher, 50); } catch {}
    showCanvas();
  }, 60);
}

// ====================================================================
//     DIRECT PROGRAMMATIC LOAD (fast path; ChemDoodle API)
// ====================================================================

function loadJsonIntoSketcher(jsonString) {
  if (typeof ChemDoodle === 'undefined' || typeof ChemDoodle.io === 'undefined') {
    throw new Error('ChemDoodle library not available.');
  }
  if (!window.sketcher) {
    throw new Error('Sketcher not ready.');
  }

  const interp = new ChemDoodle.io.JSONInterpreter();
  let contentLoaded = false;

  if (typeof interp.contentFrom === 'function') {
    try {
      const content = interp.contentFrom(jsonString);
      const mols = (content && (content.molecules || content.m || []));
      const shapes = (content && (content.shapes || content.s || []));
      const hasStuff = (Array.isArray(mols) && mols.length) || (Array.isArray(shapes) && shapes.length);
      if (hasStuff && typeof window.sketcher.loadContent === 'function') {
        window.sketcher.loadContent(mols || [], shapes || []);
        contentLoaded = true;
      }
    } catch (e) { logStep('contentFrom failed', e?.message || e); }
  }

  if (!contentLoaded && typeof interp.molFrom === 'function' && typeof window.sketcher.loadMolecule === 'function') {
    try {
      const mol = interp.molFrom(jsonString);
      if (mol) { window.sketcher.loadMolecule(mol); contentLoaded = true; }
    } catch (e) { logStep('molFrom failed', e?.message || e); }
  }

  if (!contentLoaded) throw new Error('Sketcher load failed: JSON not compatible with current ChemDoodle build.');

  window.sketcher.repaint();
  alignMoleculesLeft(window.sketcher, 50);
  showCanvas();
  logStep('sketcher repainted and aligned');
}

// ====================================================================
//     EXERCISE SELECTION HANDLER (normalize -> direct -> fallback)
// ====================================================================

async function handleExerciseSelection(exNumber) {
  if (!Number.isInteger(exNumber)) return;

  logStep('exercise selected', exNumber);

  // Optional: prevent flicker while importing
  hideCanvas();

  // Always update the exercise image/instructions
  try { showExercise(exNumber); } catch {}

  const key = String(exNumber);
  if (!EX_JSON_MAP || !EX_JSON_MAP.has(key)) {
    logStep('no JSON for exercise', { lookedFor: key, availableKeys: EX_JSON_MAP ? [...EX_JSON_MAP.keys()] : [] });
    showCanvas();
    return;
  }

  const raw = EX_JSON_MAP.get(key);
  logStep('raw CSV value (first 200)', (raw || '').slice(0, 200) + ((raw && raw.length > 200) ? '…' : ''));

  let jsonString;
  try {
    jsonString = normalizeJsonString(raw);
  } catch (e) {
    console.error('[CHEM] normalizeJsonString failed', e);
    showCanvas();
    return;
  }

  logStep('normalized JSON (first 200)', jsonString.slice(0, 200) + (jsonString.length > 200 ? '…' : ''));

  // Try direct API first
  try {
    loadJsonIntoSketcher(jsonString);
    logStep('loaded via direct API');
    return;
  } catch (e) {
    logStep('direct load failed, falling back to admin UI', e?.message || e);
  }

  // Fallback: admin importer (if present)
  try {
    await importViaAdminUI(jsonString);
  } catch (e) {
    console.error('[CHEM] admin UI import failed', e);
    showCanvas();
  }
}

// ====================================================================
//     SKETCHER READINESS GATE (same logic as the reference base.js)
// ====================================================================

function whenSketcherReady(cb) {
  if (window.sketcher) { cb(); return; }
  window.addEventListener('chemdoodle:ready', cb, { once: true });
}

// ====================================================================
//     DOMContentLoaded: keep your existing init + add autoload logic
// ====================================================================

// Initialize on DOM ready
document.addEventListener('DOMContentLoaded', function() {
  console.log('CurlyArrows initialized');

  // If exercise parameter in URL, select it (existing behavior)
  if (exerciseParam) {
    const sel = document.getElementById('exerciseSelect');
    if (sel) {
      sel.value = exerciseParam;
      const n = parseInt(exerciseParam, 10);
      if (!isNaN(n)) {
        try { showExercise(n); } catch {}
      }
    }
  }

  // Initialize clipboard functionality
  initializeClipboard();

  // NEW: once the sketcher is ready, preload CSV and load selected exercise JSON
  whenSketcherReady(async () => {
    try {
      await preloadExerciseJsonMap();
    } catch (e) {
      console.error('[CHEM] CSV load failed', e);
      // Still allow page to function without exercise JSON
      showCanvas();
      return;
    }

    const selectEl = document.getElementById('exerciseSelect');

    // Load immediately if:
    // - URL has ?exercise= (coming back from other tab), OR
    // - dropdown already has a value (restored state)
    const initial =
      (exerciseParam && String(exerciseParam).trim()) ||
      (selectEl && selectEl.value ? String(selectEl.value).trim() : '');

    if (initial) {
      const n0 = parseInt(initial, 10);
      if (!isNaN(n0)) {
        try { await handleExerciseSelection(n0); }
        catch (e) { console.error('[CHEM] initial exercise load failed', e); showCanvas(); }
      } else {
        showCanvas();
      }
    } else {
      showCanvas();
    }

    // Also wire change handler so selecting in this tab loads JSON immediately
    if (selectEl) {
      selectEl.addEventListener('change', async () => {
        const n = parseInt(selectEl.value, 10);
        if (isNaN(n)) return;
        try { await handleExerciseSelection(n); }
        catch (e) { console.error('[CHEM] handle selection', e); showCanvas(); }
      });
    }
  });
});

// ====================================================================
//     Exercise image controls (your original code, kept)
// ====================================================================

const img = document.getElementById('exerciseImage');
const select = document.getElementById('exerciseSelect');
const btn = document.getElementById('instructionsBtn');

function showExercise(n) {
  if (!img) return;

  const padded = String(n).padStart(2, '0');
  const paddedSrc = `/static/exe_img/exercise_${padded}.png`;
  const fallbackSrc = `/static/exe_img/exercise_${n}.png`;

  img.onerror = function handleMissing() {
    img.onerror = null;
    if (img.src.indexOf(paddedSrc) !== -1) {
      img.src = fallbackSrc;
    }
  };

  img.src = paddedSrc;
  img.alt = `Exercise ${n}`;

  showExerciseInstructions(n);
}

function showExerciseInstructions(n) {
  const banner = document.getElementById('exerciseInstructionBanner');
  const title = document.getElementById('exerciseInstructionTitle');
  const content = document.getElementById('exerciseInstructionContent');

  if (!banner || !title || !content) return;

  const padded = String(n).padStart(2, '0');
  const imageSrc = `/static/exe_img/exercise_${padded}.png`;

  title.textContent = `Exercise ${n} - Instructions`;
  content.innerHTML = `
    <div class="exercise-banner-content">
      <img src="${imageSrc}" alt="Exercise ${n}" class="exercise-banner-image"
           onerror="this.src='/static/exe_img/exercise_${n}.png'">
    </div>
  `;

  banner.style.display = 'block';
  window.scrollTo({ top: 0, behavior: 'smooth' });
}

if (select) {
  // NOTE: we still keep your existing image update,
  // but loading JSON is now handled inside whenSketcherReady wiring above.
  select.addEventListener('change', () => {
    const n = parseInt(select.value, 10);
    if (!isNaN(n)) showExercise(n);
  });
}

if (btn) {
  btn.addEventListener('click', () => {
    if (img) {
      img.onerror = null;
      img.src = '/static/exe_img/instructions.png';
      img.alt = 'Instructions';
    }
  });
}

// ====================================================================
//     Template/JSON functions (your original code, kept)
// ====================================================================

function createTemplateJSON() {
  if (typeof sketcher === 'undefined') {
    console.error('Sketcher not initialized');
    return;
  }

  const textTemplate = document.getElementById('templateInput');
  const base64Template = document.getElementById('base64templateInput');

  let molstring = JSON.stringify(
    new ChemDoodle.io.JSONInterpreter().contentTo(sketcher.molecules, sketcher.shapes)
  );

  if (textTemplate) textTemplate.value = molstring;
  if (base64Template) {
    base64Template.value = window.location.origin + window.location.pathname +
      '?template=' + btoa(molstring);
  }

  localStorage.setItem("templateString", molstring);
  console.log('Template JSON created');
}

function getStudentJsonCode() {
  if (typeof createTemplateJSON === 'function') {
    createTemplateJSON();
  }

  const tplEl = document.getElementById('templateInput');
  if (tplEl && tplEl.value && tplEl.value.trim()) {
    return tplEl.value.trim();
  }

  const stored = localStorage.getItem('templateString');
  if (stored) return stored;

  throw new Error('No JSON from canvas. Draw the mechanism and try again.');
}

// ====================================================================
//     Analyze button functionality (your original code, kept)
// ====================================================================

const analyzeBtn = document.getElementById('analyze-btn');
const analyzeInput = document.getElementById('analyze-input');
const analyzeOutput = document.getElementById('analyze-output');
const exerciseSelect = document.getElementById('exerciseSelect');

if (analyzeBtn) {
  analyzeBtn.addEventListener('click', async () => {
    try {
      if (analyzeOutput) analyzeOutput.value = '';

      if (!exerciseSelect || !exerciseSelect.value) {
        if (analyzeOutput) analyzeOutput.value = '❌ Please select an exercise first!';
        return;
      }

      const exercise = parseInt(exerciseSelect.value, 10);
      const student_reasoning = (analyzeInput ? analyzeInput.value : '').trim();
      const student_json_code = getStudentJsonCode();

      if (analyzeOutput) {
        analyzeOutput.value = '⏳ Analyzing your mechanism...\n\n' +
          'Exercise: ' + exercise + '\n' +
          'Reasoning length: ' + (student_reasoning.length) + ' characters\n' +
          'Mechanism data: ' + (student_json_code.length) + ' characters';
      }
      analyzeBtn.disabled = true;
      analyzeBtn.innerHTML = '<span class="spinner-border spinner-border-sm"></span> Analyzing...';

      const res = await fetch('/api/analyze', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        credentials: 'include',
        body: JSON.stringify({
          exercise,
          student_json_code,
          student_reasoning: student_reasoning || null
        })
      });

      if (!res.ok) {
        const errorText = await res.text();
        throw new Error(errorText || 'Server error');
      }

      const data = await res.json();

      if (analyzeOutput) {
        analyzeOutput.value =
          '✅ Analysis Complete\n\n' +
          (typeof data === 'string'
            ? data
            : data.feedback || data.result || JSON.stringify(data, null, 2));

        if (data.submission_id) {
          analyzeOutput.value += '\n\n--- Submission ID: ' + data.submission_id + ' ---';
        }
      }

    } catch (e) {
      console.error('Analysis error:', e);
      if (analyzeOutput) analyzeOutput.value = '❌ Error: ' + (e?.message || e);
    } finally {
      if (analyzeBtn) {
        analyzeBtn.disabled = false;
        analyzeBtn.innerHTML = '<i class="bi bi-cpu"></i> <strong>ANALYZE</strong>';
      }
    }
  });
}

// ====================================================================
//     Clipboard functionality (your original code, kept)
// ====================================================================

function initializeClipboard() {
  const copyBtn = document.getElementById('copyToClipboard');
  const clearBtn = document.getElementById('clearClipboard');
  const saveBtn = document.getElementById('saveClipboard');

  if (clearBtn) {
    clearBtn.addEventListener('click', () => {
      if (typeof sketcher !== 'undefined') {
        window.location.reload();
      }
    });
  }

  if (copyBtn) {
    copyBtn.addEventListener('click', async () => {
      try {
        if (typeof sketcher === 'undefined' || typeof ChemDoodle === 'undefined') {
          alert('Sketcher not initialized');
          return;
        }

        const chemimage = ChemDoodle.io.png.string(sketcher);
        const response = await fetch(chemimage);
        const blob = await response.blob();

        await navigator.clipboard.write([
          new ClipboardItem({ [blob.type]: blob })
        ]);

        createTemplateJSON();
        alert('Copied to clipboard!');

      } catch (error) {
        console.error('Copy error:', error);
        alert('Failed to copy: ' + error.message);
      }
    });
  }

  if (saveBtn) {
    saveBtn.addEventListener('click', () => {
      if (typeof sketcher === 'undefined' || typeof ChemDoodle === 'undefined') {
        alert('Sketcher not initialized');
        return;
      }

      const chemimage = ChemDoodle.io.png.string(sketcher);
      const now = new Date();
      const date = now.toISOString().split("T")[0];
      const time = now.toISOString().split("T")[1];
      const filename = date.replace(/-/g, "") + time.replace(/:/g, "").split(".")[0] + ".png";

      const hiddenElement = document.createElement("a");
      hiddenElement.href = chemimage;
      hiddenElement.target = "_blank";
      hiddenElement.download = filename;
      hiddenElement.click();

      createTemplateJSON();
    });
  }
}

// Timeline/History functionality (optional enhancement)
const openTimelineBtn = document.querySelector('.openTimeline');
const closeTimelineBtn = document.querySelector('.closeTimeline');
const timeline = document.querySelector('.timeline');

if (openTimelineBtn && timeline) {
  openTimelineBtn.addEventListener('click', () => {
    timeline.classList.add('open');
  });
}

if (closeTimelineBtn && timeline) {
  closeTimelineBtn.addEventListener('click', () => {
    timeline.classList.remove('open');
  });
}

console.log('Base.js loaded successfully (UPDATED)');