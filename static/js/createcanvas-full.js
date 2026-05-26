ChemDoodle.ELEMENT["H"].jmolColor = "black";
ChemDoodle.ELEMENT["S"].jmolColor = "#B9A130";

// Calculate responsive canvas size
// Get the actual container element to match its width precisely
var container = document.getElementById("sketcherContainer");
var containerWidth = container ? container.offsetWidth - 30 : window.innerWidth - 30; // Subtract padding (15px * 2)
var canvasWidth = containerWidth; // Use full container width
var canvasHeight = Math.min(400, window.innerHeight * 0.6);

var sketcher = new ChemDoodle.SketcherCanvas("sketcher", canvasWidth, canvasHeight, {
  useServices: false,
  resizable: false,
});
sketcher.styles.atoms_displayTerminalCarbonLabels_2D = true;
sketcher.styles.atoms_useJMOLColors = true;
sketcher.styles.bonds_clearOverlaps_2D = true;
sketcher.styles.shapes_color = "#c10000";

// *****ButtonDisplay*****
// console.log('url search: ' + window.location.search)

function setButtons(){
  const parameters = window.location.search;
  const parametersArray = parameters.split(/[?&]/)
  const adminTest = (element) => element === 'admin=true';
  let resultAdminTest = parametersArray.some(adminTest)
  console.log(`Setbuttons - Admintest: ${resultAdminTest}`)
  let invisibleBtn = [
    // Add the id elements of the as a string
    // id elements can be listed as json format example:
    // 'idname1',
    // 'idname2',
    // etc.
    // 'sketcher_button_open',
    'sketcher_button_save'
  ];
  if(resultAdminTest == false){
    invisibleBtn.push('sketcher_button_template_label')
  }
  // Pushing buttons to Array
  invisibleBtn.forEach(button => {
    document.getElementById(button).style.display = "none";
    document.getElementById(button).style.visibility = "hidden";
  });
}
setButtons();
// Repaint canvas
sketcher.repaint();

