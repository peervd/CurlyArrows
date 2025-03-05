from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator,rdFMCS, BRICS
from rdkit.DataStructs import TanimotoSimilarity



def compute_substructure_score(mol1, mol2, timeout=1):
    """
    Compute a substructure similarity score based on the maximum common substructure (MCS).
    The score is defined as:
         score = (number of atoms in MCS) / (min(number of atoms in mol1, mol2))
    If no MCS is found, the score is 0.
    """
    res = rdFMCS.FindMCS([mol1, mol2], timeout=timeout)
    if res.canceled or res.numAtoms == 0:
        return 0.0
    return res.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())

def combined_similarity(mol1, mol2, fp1, fp2):
    """
    Combine the Tanimoto similarity of Morgan fingerprints and an MCS-based substructure score.
    Here we take the maximum of the two scores.
    """
    t_sim = TanimotoSimilarity(fp1, fp2)
    s_sim = compute_substructure_score(mol1, mol2)
    return max(t_sim, s_sim)

def get_candidate_fragments(smi, min_atoms=4):
    """
    For a given molecule SMILES, return a set of candidate fragments.
    Always include the parent canonical SMILES.
    If the molecule has at least `min_atoms` heavy atoms,
    add its BRICS fragments (if any).
    """
    candidates = set()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return candidates
    parent = Chem.MolToSmiles(mol, canonical=True)
    candidates.add(parent)
    
    # Only fragment if the molecule is sufficiently large.
    if mol.GetNumHeavyAtoms() >= min_atoms:
        try:
            frags = BRICS.BRICSDecompose(mol)
            # BRICSDecompose returns a set of fragment SMILES.
            candidates.update(frags)
        except Exception as e:
            pass
    return candidates

def canonicalize_reaction(smiles_reaction, sim_threshold=0.5):
    """
    Convert a reaction SMILES into a canonical form by removing spectator molecules.
    
    This function first fragments each reactant and product (using BRICS) to generate candidate
    substructures. Then, for each parent molecule, if any of its candidate fragments finds a matching
    candidate on the opposite side (with combined similarity >= sim_threshold), the parent is retained.
    
    Molecules that do not find any matching partner are considered spectators and removed.
    
    Parameters:
      smiles_reaction (str): Reaction SMILES with reactants and products separated by '>>'
                             and individual molecules separated by '.'.
      sim_threshold (float): Minimum similarity (0–1) required for a candidate match.
    
    Returns:
      str: Canonical reaction SMILES including only the core transformation.
    """
    try:
        reactants_str, products_str = smiles_reaction.split(">>")
    except Exception as e:
        raise ValueError("Reaction SMILES must contain '>>' separating reactants and products.")
    
    reactant_smiles = reactants_str.split(".")
    product_smiles = products_str.split(".")
    
    # For each parent molecule, generate candidate fragments.
    reactant_candidates = {}  # key: parent SMILES, value: set of candidate SMILES
    product_candidates = {}
    
    for smi in reactant_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            reactant_candidates[smi] = cand
    for smi in product_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            product_candidates[smi] = cand
    
    # Create a fingerprint generator.
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    # For each candidate fragment, precompute its RDKit molecule and fingerprint.
    def get_fp(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None, None
        fp = generator.GetFingerprint(mol)
        return mol, fp
    
    reactant_frag_data = {}  # key: candidate SMILES, value: (mol, fp)
    for parent, cand_set in reactant_candidates.items():
        for cand in cand_set:
            mol, fp = get_fp(cand)
            if mol is not None:
                reactant_frag_data[cand] = (mol, fp)
                
    product_frag_data = {}
    for parent, cand_set in product_candidates.items():
        for cand in cand_set:
            mol, fp = get_fp(cand)
            if mol is not None:
                product_frag_data[cand] = (mol, fp)
    
    # Determine which parent molecules (reactant and product) are "involved".
    involved_reactants = set()
    involved_products = set()
    
    # For each reactant parent, if any of its candidate fragments finds a match in any product candidate, mark it as involved.
    for r_parent, r_cands in reactant_candidates.items():
        for r_cand in r_cands:
            if r_cand not in reactant_frag_data:
                continue
            r_mol, r_fp = reactant_frag_data[r_cand]
            for p_parent, p_cands in product_candidates.items():
                for p_cand in p_cands:
                    if p_cand not in product_frag_data:
                        continue
                    p_mol, p_fp = product_frag_data[p_cand]
                    sim = combined_similarity(r_mol, p_mol, r_fp, p_fp)
                    if sim >= sim_threshold:
                        involved_reactants.add(r_parent)
                        involved_products.add(p_parent)
                        # Once one match is found for this reactant parent, no need to check further.
                        break
                else:
                    continue
                break
    
    # Also check from the product side (for safety)
    for p_parent, p_cands in product_candidates.items():
        for p_cand in p_cands:
            if p_cand not in product_frag_data:
                continue
            p_mol, p_fp = product_frag_data[p_cand]
            for r_parent, r_cands in reactant_candidates.items():
                for r_cand in r_cands:
                    if r_cand not in reactant_frag_data:
                        continue
                    r_mol, r_fp = reactant_frag_data[r_cand]
                    sim = combined_similarity(r_mol, p_mol, r_fp, p_fp)
                    if sim >= sim_threshold:
                        involved_reactants.add(r_parent)
                        involved_products.add(p_parent)
                        break
                else:
                    continue
                break
    
    # Build the canonical reaction SMILES from the parent molecules that are involved.
    core_reactants_smiles = ".".join(sorted(involved_reactants))
    core_products_smiles = ".".join(sorted(involved_products))
    return f"{core_reactants_smiles}>>{core_products_smiles}"


def find_longest_matching_sequences(model_canonical, student_canonical):
    """Find the longest non-overlapping matching subsequences between the model and student's sequences, including shift."""
    model_len = len(model_canonical)
    student_len = len(student_canonical)
    
    matches = []
    
    for i in range(model_len):
        for j in range(i + 1, model_len + 1):
            sub_seq = model_canonical[i:j]
            for k in range(student_len - len(sub_seq) + 1):
                if student_canonical[k:k+len(sub_seq)] == sub_seq:
                    shift = k - i  # Calculate shift
                    matches.append((i, j - 1, shift))
    
    # Filter out overlapping sequences, keeping only the longest ones
    non_overlapping_matches = []
    matches.sort(key=lambda x: x[1] - x[0], reverse=True)  # Sort by length, longest first
    
    for start, end, shift in matches:
        if not any(s <= start <= e or s <= end <= e for s, e, _ in non_overlapping_matches):
            non_overlapping_matches.append((start, end, shift))
    
    return sorted(non_overlapping_matches)

def compare_reactions(model_list, student_list):
    """Compare a model reaction sequence with a student's reaction sequence.
    
    The "individual_steps" key in the output dictionary represents a list of 
    Boolean values. Each entry in this list corresponds to a reaction in the 
    student's list and indicates whether that specific reaction appears in the 
    model list after canonicalization (i.e., after removing spectator 
    molecules)
    
    The "matching_sequences" output contains a list of tuples, where each tuple 
    (i,j) represents a subsequence of the model list that is also present in 
    the student's list. Specifically, it means that the transformations from 
    index i to j (inclusive) in the model list appear somewhere in the student 
    list in the correct order, but possibly shifted in position.

    For example, if "matching_sequences": [(5, 10, 0)], this means that 
    transformations from step 5 to step 10 in the model sequence appear 
    somewhere in the student's list in the correct order. Reactants are always
    indexed as 0.
    """
    model_canonical = [canonicalize_reaction(rx) for rx in model_list]
    student_canonical = [canonicalize_reaction(rx) for rx in student_list]
    individual_comparisons = [student_rx in model_canonical for student_rx in student_canonical]
    matching_subsequences = find_longest_matching_sequences(model_canonical, student_canonical)
    
    return {
        "individual_steps": individual_comparisons,
        "matching_sequences": matching_subsequences
    }
