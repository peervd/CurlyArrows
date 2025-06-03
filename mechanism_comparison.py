import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from exersices_ORC import load_acid_base
from helper_functions import count_trues, _split_reaction,all_lists_empty
from molecular_structures import sanitize_hypervalent_smiles, create_mol_safely, combined_similarity, get_candidate_fragments, compare_molecules

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


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
        reactant_smiles,product_smiles = _split_reaction(smiles_reaction)
    except Exception as e:
        raise ValueError("Reaction SMILES must contain '>>' separating reactants and products.")

    # For each parent molecule, generate candidate fragments.
    reactant_candidates = {}  # key: parent SMILES, value: set of candidate SMILES
    product_candidates = {}
    unsanitized_r_mol = []
    unsanitized_p_mol = []
    for smi in reactant_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            reactant_candidates[smi] = cand
        else:
            unsanitized_r_mol.append(smi)
    for smi in product_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            product_candidates[smi] = cand
        else:
            unsanitized_p_mol.append(smi)

    # Add the unsanitized molecules back as their own candidates
    # This ensures they are still considered in the reaction
    for smi in unsanitized_r_mol:
        reactant_candidates[smi] = {smi}
    for smi in unsanitized_p_mol:
        product_candidates[smi] = {smi}

    # Create a fingerprint generator.
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    # For each candidate fragment, precompute its RDKit molecule and fingerprint.
    def get_fp(smi):
        mol = create_mol_safely(smi)
        if mol is None:
            return None, None
        try:
            fp = generator.GetFingerprint(mol)
            return mol, fp
        except:
            return None, None
    
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
    
    # String-based fallback for molecules that couldn't be properly converted to RDKit mols
    # Get raw smiles that haven't been classified as involved yet
    reactant_raw_smiles = {smi for smi in reactant_smiles if smi not in involved_reactants}
    product_raw_smiles = {smi for smi in product_smiles if smi not in involved_products}
    
    # This is a simple direct comparison of SMILES, which works for exact matches
    for r_smi in reactant_raw_smiles:
        sanitized_r = sanitize_hypervalent_smiles(r_smi)
        for p_smi in product_raw_smiles:
            sanitized_p = sanitize_hypervalent_smiles(p_smi)
            # Simple string comparison as fallback
            if sanitized_r == sanitized_p:
                involved_reactants.add(r_smi)
                involved_products.add(p_smi)
    reactants = []
    for x in involved_reactants:
        try:
            c = Chem.MolFromSmiles(x)
            c = Chem.MolToSmiles(c)
            reactants.append(c)
        except:
            reactants.append(x)
    products = []
    for x in involved_products:
        try:
            c = Chem.MolFromSmiles(x)
            c = Chem.MolToSmiles(c)
            products.append(c)
        except:
            products.append(x)
     
    # Build the canonical reaction SMILES from the parent molecules that are involved.
    core_reactants_smiles = ".".join(sorted(reactants))
    core_products_smiles = ".".join(sorted(products))
    return f"{core_reactants_smiles}>>{core_products_smiles}"


def find_longest_matching_sequences(steps_dict):
    """Find the longest non-overlapping matching subsequences between the model and student's sequences, including shift."""
    steps = list(steps_dict.keys())
    missing_steps = []
    new_dict = {}
    for x in steps:
        if len(steps_dict[x]) > 1:
            length = len(steps_dict[x])
            if steps_dict[x][-1][0][1] == True:   # check if last step is True (for example product, or further reaction intermediate)
                new_dict[x] = ([True,True],steps_dict[x][0][1])
                miss_arr = np.arange(x+1,(x+length),1)
                missing_steps.extend(miss_arr.tolist())
            elif steps_dict[x][-1][0][1] == False:
                new_dict[x] = ([True,False],steps_dict[x][0][1])
                miss_arr = np.arange(x+2,(x+length),1)
                missing_steps.extend(miss_arr.tolist())
        elif steps_dict[x]:
            new_dict[x] = steps_dict[x][0]
        elif not steps_dict[x]:
            steps.remove(x)
            del steps_dict[x]
    if all(new_dict[s][1] == 'not_present' for s in list(new_dict.keys())):
        non_overlapping_matches = []
        return non_overlapping_matches
    matches = []

    for i,step in enumerate(steps):
        if new_dict[step][1] == ('not_present'):
            matches.append(('not_present',0))
        elif new_dict[step][1] == ('no_chemical_transformations'):
            pass
        else:
            shift = step - new_dict[step][1]
            if new_dict[step][0] == [True,True]:
                matches.append((step,shift))
            elif new_dict[step][0] == [True,False]:
                matches.append((step,shift))
            elif new_dict[step][0] == [False,True]:
                matches.append((step+1,shift))

    # Initialize tracking variables
    current_start = None
    current_end = None
    current_shift = None
    non_overlapping_matches = []


    for i, (step, shift) in enumerate(matches):
        # Helper function to check if values are consecutive integers
        def is_consecutive(a, b):
            if isinstance(a, str) or isinstance(b, str):
                return False
            return a == b + 1
    
        if current_start is None:
            # Start a new sequence
            current_start = step
            current_end = step
            current_shift = shift
        # Check if current step continues the sequence (both are integers and consecutive)
        elif step != 'not_present' and current_end != 'not_present' and is_consecutive(step, current_end):
            current_end = step 
        else:
            # Either we have a 'not_present' or values aren't consecutive
            if current_start != 'not_present':
                new_sequence = (current_start, current_end, current_shift)
                non_overlapping_matches.append(new_sequence)
            
            # Start a new sequence if step is not 'not_present'
            if step != 'not_present':
                current_start = step
                current_end = step
                current_shift = shift
            else:
                current_start = None
                current_end = None
                current_shift = None

    if current_start != None and current_end != None and current_start != 'not_present' and current_end != 'not_present':
        only_sequence = (current_start, current_end, current_shift)
        non_overlapping_matches.append(only_sequence)

    return sorted(non_overlapping_matches),missing_steps

def compare_steps(model_canonical, student_canonical, acid_base):
    """
    Compare two lists of chemical reactions and determine how student reactions match with model reactions.
    
    Parameters:
    - model_canonical: List of model reactions in SMILES format
    - student_canonical: List of student reactions in SMILES format
    - acid_base: List of acid/base pairs that can be substituted
    
    Returns:
    - Dictionary with indices of student reactions as keys and matching information as values
    """
    # Get the lengths of both lists
    st_len = len(student_canonical)
    m_len = len(model_canonical)
    
    # Initialize the results dictionaries
    individual_comparisons = {}
    molecular_structure = {}
    product = False
    # Parse all reactions once to avoid repeated parsing
    model_reactions = []
    for m_rx in model_canonical:
        parts = m_rx.split('>>')
        m_r = parts[0].split('.')
        m_p = parts[1].split('.')
        model_reactions.append((m_r, m_p))
    
    student_reactions = []
    for st_rx in student_canonical:
        parts = st_rx.split('>>')
        st_r = parts[0].split('.')
        st_p = parts[1].split('.')
        student_reactions.append((st_r, st_p))

    # Determine which approach to use based on which list is longer
    if st_len >= m_len:
        # Case: Student has more or equal reactions compared to model
        individual_comparisons = {key: [] for key in range(st_len)}
        molecular_structure = {key: None for key in range(st_len)}
        
        for st_idx, (st_r, st_p) in enumerate(student_reactions):
            reaction_found = False
            if st_r == [''] and st_p == ['']:
                individual_comparisons[st_idx] = [([False, False], 'no_chemical_transformations')]
                continue
            # Compare this student reaction with each model reaction
            for m_idx, (m_r, m_p) in enumerate(model_reactions):
                # Check reactants and products
                r_match_info = check_molecules_match(st_r, m_r, acid_base)
                p_match_info = check_molecules_match(st_p, m_p, acid_base)
                
                bool_r = r_match_info['all_match']
                bool_p = p_match_info['all_match']
                if m_len-1 == m_idx and bool_p == True:
                    product = True
                    
                # Determine the overall match status
                if bool_r and bool_p:
                    individual_comparisons[st_idx] = [([True, True], m_idx)]
                    
                    
                    reaction_found = True
                    break
                elif bool_r:
                    individual_comparisons[st_idx].append(([True, False], m_idx))
                    if True in p_match_info['individual_matches'] and False in p_match_info['individual_matches'] and p_match_info['unique_molecules'] and p_match_info['unique_from_model']:
                        differences = analyze_molecular_differences(
                            p_match_info['unique_molecules'], 
                            p_match_info['unique_from_model'])
                        if all(x == None for x in differences):
                            pass
                        else:
                            molecular_structure[st_idx+1] = differences
                        
                    reaction_found = True
                    continue
                elif bool_p:
                    individual_comparisons[st_idx].append(([False, True], m_idx))

                    reaction_found = True
                    # If there's a mismatch, analyze the molecular differences
                    
                    if True in r_match_info['individual_matches'] and False in r_match_info['individual_matches'] and r_match_info['unique_molecules'] and r_match_info['unique_from_model']:
                        if molecular_structure[st_idx] != None:
                            pass
                        else:
                            differences = analyze_molecular_differences(
                                r_match_info['unique_molecules'], 
                                r_match_info['unique_from_model'])
                            if all(x == None for x in differences):
                                pass    
                            else:
                                molecular_structure[st_idx] = differences
                        
                    continue
            
            # If no match found for this student reaction
            if not reaction_found:
                individual_comparisons[st_idx].append(([False, False], 'not_present'))
                
    
    else:
        # Case: Model has more reactions than student
        individual_comparisons = {key: [] for key in range(st_len)}
        molecular_structure = {key: None for key in range(st_len)}

        for m_idx, (m_r, m_p) in enumerate(model_reactions):
            reaction_found = False

            for st_idx, (st_r, st_p) in enumerate(student_reactions):
                # Skip if we already found a perfect match for this student reaction
                if st_idx in individual_comparisons and individual_comparisons[st_idx] and individual_comparisons[st_idx][0][0] == [True, True]:
                    continue
                
                # Check reactants and products
                r_match_info = check_molecules_match(st_r, m_r, acid_base)
                p_match_info = check_molecules_match(st_p, m_p, acid_base)
                bool_r = r_match_info['all_match']
                bool_p = p_match_info['all_match']
                if m_len-1 == m_idx and bool_p == True:
                    product = True
                    
                # Determine the overall match status
                if bool_r and bool_p:
                    individual_comparisons[st_idx] = [([True, True], m_idx)]
                    
                    reaction_found = True
                    break
                elif bool_r and (st_idx not in individual_comparisons or individual_comparisons[st_idx] != [([True, True], m_idx)]):
                    individual_comparisons[st_idx].append(([True, False], m_idx))
                    if True in p_match_info['individual_matches'] and False in p_match_info['individual_matches'] and p_match_info['unique_molecules'] and p_match_info['unique_from_model']:
                        differences = analyze_molecular_differences(
                            p_match_info['unique_molecules'], 
                            p_match_info['unique_from_model'])
                        if all(x == None for x in differences):
                            pass    
                        else:
                            molecular_structure[st_idx+1] = differences
                    
                    reaction_found = True                   
                    continue
                elif bool_p and (st_idx not in individual_comparisons or individual_comparisons[st_idx] != [([True, True], m_idx)]):
                    individual_comparisons[st_idx].append(([False, True], m_idx))
                    
                    reaction_found = True
                    # If there's a mismatch, analyze the molecular differences
                    if True in r_match_info['individual_matches'] and False in r_match_info['individual_matches'] and r_match_info['unique_molecules'] and r_match_info['unique_from_model']:
                        if st_idx in list(molecular_structure.keys()):
                            pass
                        else:
                            differences = analyze_molecular_differences(
                                r_match_info['unique_molecules'], 
                                r_match_info['unique_from_model'])
                            if all(x == None for x in differences):
                                pass    
                            else:
                                molecular_structure[st_idx] = differences 
                    continue
            # If no match found for this model reaction
            if not reaction_found:
                if all_lists_empty(individual_comparisons):
                    individual_comparisons[0].append(([False, False], "not_present"))
                else:
                    ind_add = max([k for k, v in individual_comparisons.items() if v]) # find key in dict with a not empty list
                    individual_comparisons[ind_add].append(([False, False], "not_present"))
                
    # Clean up the molecular structure dictionary
    molecular_structure = {k: v for k, v in molecular_structure.items() if v is not None}
    
    # Add general information
    general = {
        "model_steps": m_len,
        "student_steps": st_len,
    }
    ### remove empty lists as values from individual_comparisons
    individual_comparisons = {k: v for k, v in individual_comparisons.items() if v}
    
    return individual_comparisons, molecular_structure, general, {"product":product}


def check_molecules_match(molecules1, molecules2, acid_base):
    """
    Check how molecules from list1 match with molecules from list2, considering acid/base substitutions.
    
    Returns a dictionary with match information.
    """
    individual_matches = []
    matching_molecules = []
    unique_molecules = []
    
    for mol1 in molecules1:
        if mol1 in molecules2:
            individual_matches.append(True)
            matching_molecules.append(mol1)
        elif any(mol1 in acid[0] for acid in acid_base) and any(acid[0] in molecules2 for acid in acid_base):
            individual_matches.append(True)
            matching_molecules.append(mol1)
        elif any(mol1 in base[1] for base in acid_base) and any(base[1] in molecules2 for base in acid_base):
            individual_matches.append(True)
            matching_molecules.append(mol1)
        else:
            individual_matches.append(False)
            unique_molecules.append(mol1)
    
    unique_from_model = set(molecules2) - set(matching_molecules)
    
    return {
        'individual_matches': individual_matches,
        'all_match': all(individual_matches),
        'matching_molecules': matching_molecules,
        'unique_molecules': unique_molecules,
        'unique_from_model': unique_from_model
    }


def analyze_molecular_differences(unique_molecules, unique_from_model):
    """
    Analyze differences between molecules that don't match directly.
    """
    if not unique_from_model:
        # If there's no matching molecule in the model, report it
        return ['molecule_not_present_in_model', sum(1 for atom in Chem.MolFromSmiles(unique_molecules[0]).GetAtoms() if atom.GetSymbol() == 'C')]
    
    differences = []
    for mol in unique_molecules:
        comparison = []
        for model_mol in unique_from_model:
            comparison.append(compare_molecules(mol, model_mol))
        
        # Find the best match based on priority
        pri = 5
        diff = None
        for z in comparison:
            if z[0].endswith('_1') and pri > 0:
                diff = [z[0].strip('_1'), z[1]]
                pri = 0
            elif z[0].endswith('_2') and pri > 1:
                diff = [z[0].strip('_2'), z[1]]
                pri = 1
            elif z[0].endswith('_3') and pri > 2:
                diff = [z[0].strip('_3'), z[1]]
                pri = 2
            elif z[0].endswith('_4') and pri > 3:
                diff = [z[0].strip('_4'), z[1]]
                pri = 3
            elif z[0].endswith('_0') and pri > 4:
                pri = 4
        
        differences.append(diff)
    
    return differences



def compare_reactions(model_list, student_list):
    """
    Compare a model reaction sequence with a student's reaction sequence.
    
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
    student_canonical = [canonicalize_reaction(rx) for rx in student_list]
    acid_base = load_acid_base('list_acid_base.csv')
    index_comp = False

    if isinstance(model_list, list) and all(isinstance(item, list) for item in model_list):
        lst_comparisons = []
        mol_found = []
        for m_lst in model_list:
            model_canonical = [canonicalize_reaction(rx) for rx in m_lst]
            ind_comp = compare_steps(model_canonical,student_canonical,acid_base)
            lst_comparisons.append(ind_comp)
                # Count keys in the second dictionary of each tuple
        key_counts = [len(comparison[1]) for comparison in lst_comparisons]
        
        # Check if there's a tie for maximum number of keys
        max_count = max(key_counts)
        has_tie = key_counts.count(max_count) > 1
        
        if has_tie:
            tied_indices = [i for i, count in enumerate(key_counts) if count == max_count]
            
            # Select the first one with max count if you need to choose one
            index_comp = key_counts.index(max_count)
            selected_comparison = lst_comparisons[index_comp]
            individual_comparisons = lst_comparisons[index_comp]
        else:
            index_comp, individual_comparisons = max(enumerate(lst_comparisons), key=lambda x: count_trues(x[1]))

    else:   
        model_canonical = [canonicalize_reaction(rx) for rx in model_list]
        individual_comparisons = compare_steps(model_canonical, student_canonical,acid_base)

    matching_subsequences = find_longest_matching_sequences(individual_comparisons[0])

    if not matching_subsequences:
        return {"individual_steps": individual_comparisons[0],"matching_sequences": [],
                "missing_steps":[],"molecular_structures":individual_comparisons[1],"index_resonance": index_comp, "product": individual_comparisons[3]["product"]}
    else:
        return {"individual_steps": individual_comparisons[0],"matching_sequences": matching_subsequences[0],
                "missing_steps":matching_subsequences[1],"molecular_structures":individual_comparisons[1],"index_resonance": index_comp, "product": individual_comparisons[3]["product"]}