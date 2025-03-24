from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator,rdFMCS, BRICS
from rdkit.DataStructs import TanimotoSimilarity
from exersices_ORC import load_acid_base



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


def find_longest_matching_sequences(steps_dict):
    """Find the longest non-overlapping matching subsequences between the model and student's sequences, including shift."""
    steps = steps_dict.keys()
    if all(steps_dict[s][1] == 'not_present' for s in steps):
        non_overlapping_matches = []
        return non_overlapping_matches
    matches = []

    for i,step in enumerate(steps):
        
        shift = step - steps_dict[step][1]
        if steps_dict[step][0] == [True,True]:
            matches.append((step,shift))
        elif steps_dict[step][0] == [True,False]:
            matches.append((step,shift))
        elif steps_dict[step][0] == [False,True]:
            matches.append((step+1,shift))
    non_overlapping_matches = []
    
    # Initialize tracking variables
    current_start = None
    current_end = None
    current_shift = None
    non_overlapping_matches = []

    for i, (step, shift) in enumerate(matches):

        if current_start is None:
            # Start a new sequence
            current_start = step
            current_end = step
            current_shift = shift
        elif step == current_end + 1:
            current_end = step 
        else:
            # Store the completed sequence if it's not already present
            new_sequence = (current_start, current_end, current_shift)
            non_overlapping_matches.append(new_sequence)
            
            # Start a new sequence
            current_start = step 
            current_end = step + shift
            current_shift = shift
    
    print(non_overlapping_matches)
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
    acid_base = load_acid_base('list_acid_base.csv')

    #individual_comparisons = [student_rx in model_canonical for student_rx in student_canonical]
    individual_comparisons = {}
    for st_numb,st_rx in enumerate(student_canonical):
        st_rx = st_rx.split('>>')
        st_r = st_rx[0].split('.')
        st_p = st_rx[1].split('.')
        reaction_found = False
        for m_numb,m_rx in enumerate(model_canonical):
            m_rx = m_rx.split('>>')
            m_r = m_rx[0].split('.')
            m_p = m_rx[1].split('.')
            bool_r = []
            r = []
            bool_p = []
            p = []
            if any(x in m_r for x in st_r):
                for x in st_r:
                    if x in m_r:
                        bool_r.append(True)
                        reactant = True
                        r.append(x)
                    
                    elif [x in acid[0] for acid in acid_base] and any(acid[0] in m_r for acid in acid_base):
                        bool_r.append(True)
                        reactant = True
                        r.append(x)
                    
                    elif [x in base[1] for base in acid_base] and any(base[1] in m_r for base in acid_base):
                        bool_r.append(True)
                        r.append(x)
                        reactant = True
                    
            else: 
                bool_r.append(False)

            if any(x in m_p for x in st_p):
                for x in st_p:
                    if x in m_p:
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
                    elif [x in acid[0] for acid in acid_base] and any(acid[0] in m_p for acid in acid_base):
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
                    elif [x in base[1] for base in acid_base] and any(base[1] in m_p for base in acid_base):
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
            else: 
                bool_p.append(False)

            if bool_r and bool_p and all(x == True for x in bool_r) and all(x == True for x in bool_p):
                individual_comparisons[st_numb] = ([True,True],m_numb)
                reaction_found = True
                break

            elif bool_r and bool_p and all(x == True for x in bool_r) and all(x == False for x in bool_p):
                max_length = max([len(s) for s in m_r])
                if max_length == max([len(s) for s in r]):
                    individual_comparisons[st_numb] = ([True,False],m_numb)
                    reaction_found = True
                    break

            elif bool_r and bool_p and all(x == False for x in bool_r) and all(x == True for x in bool_p):
                max_length = max([len(s) for s in m_p])
                if max_length == max([len(s) for s in p]):
                    individual_comparisons[st_numb] = ([False,True],m_numb)
                    reaction_found = True
                    break
        if reaction_found == False:
            individual_comparisons[st_numb] = ([False,False],'not_present')
    print(individual_comparisons)
    matching_subsequences = find_longest_matching_sequences(individual_comparisons)

    return {"individual_steps": individual_comparisons,"matching_sequences": matching_subsequences}