import numpy as np
import copy
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from .helper_functions import count_trues, _split_reaction,all_lists_empty, is_consecutive
from .molecular_structures import calculate_substructre_score,safe_mol_from_smiles,safe_smiles, combined_similarity, get_candidate_fragments, compare_molecules, get_molecule_name, molecule_core, composition_counter_diff, optimize_reaction, mols_equal, has_metal, zero_charge, generate_valid_combinations, is_valid_combination, calculate_mass_balance,calculate_charge_balance, score_reaction, check_subfragment_match

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def canonicalize_reaction_list(list_reactions, sim_threshold = 0.5,mass_limit = 13):
    """
    Convert a reaction SMILES into a canonical form by removing spectator molecules for each step.
    
    This function first fragments each reactant and product for a step (using BRICS) to generate candidate
    substructures. Then, for each parent molecule, if any of its candidate fragments finds a matching
    candidate on the opposite side (with combined similarity >= sim_threshold), the parent is retained.
    Molecules with the highest match are prioritized. 
    
    Molecules that do not find any matching partner are considered spectators and removed and reported as rogue.
    
    Parameters:
      smiles_reaction (str): Reaction SMILES with reactants and products separated by '>>'
                             and individual molecules separated by '.'.
      sim_threshold (float): Minimum similarity (0–1) required for a candidate match.
    
    Returns:
      lst of str: Canonical reaction SMILES including only the core transformation.
      dct of str: True spectator molecules
    """
    
    list_can_r = []
    dct_rogue = {}
    rogue_molecules = False
    for i,reaction in enumerate(list_reactions):
        try:
            unsanitized_reactant_smiles,unsanitized_product_smiles = _split_reaction(reaction)
        except Exception as e:
            raise ValueError("Reaction SMILES must contain '>>' separating reactants and products.")
        spectators = {'r': [],'p': []}
        reactants = []
        products = []
        
        reactant_smiles = [safe_smiles(x) for x in unsanitized_reactant_smiles]
        product_smiles = [safe_smiles(x) for x in unsanitized_product_smiles]
        
        ### create dict with the core of each molecule, sanitized without H and charges        
        # contains list of core string, mol obj, charges
        mol_r = {r: safe_mol_from_smiles(r) for r in reactant_smiles}
        mol_p = {p: safe_mol_from_smiles(p) for p in product_smiles}
        
        ### check for substructure match between molecules before and after the reaction arrow
        r_match = defaultdict(list)
        p_match = defaultdict(list)
        for r in mol_r.keys():
            for p in mol_p.keys():
                mol_found = False
                if mols_equal(Chem.MolFromSmiles(r),Chem.MolFromSmiles(p)):
                    spectators['r'].append(r)
                    spectators['p'].append(p)
                    break
                elif r in spectators['r'] or p in spectators['p']:
                    break
                r0 = zero_charge(mol_r[r]) 
                p0 = zero_charge(mol_p[p])
                if r0.HasSubstructMatch(p0):
                    r_match[r].append(p)
                    mol_found = True
                elif has_metal(r0) and composition_counter_diff(p0,r0):
                    p_match[p].append(r)
                    mol_found = True
                if p0.HasSubstructMatch(r0):
                    p_match[p].append(r)
                    mol_found = True
                elif has_metal(p0) and composition_counter_diff(r0,p0):
                    p_match[p].append(r)
                    mol_found = True
                if mol_found == False and check_subfragment_match(r,p) == True:
                    r_match[r].append(p)
            if not r_match[r]:
                r_match[r] = []

        for r in spectators['r']:
            for chem in r_match.keys():
                if r == chem:
                    try:
                        r_match.pop(r)
                    except:
                        continue
                    break
                elif r in r_match[chem]:
                    try:
                        r_match.pop(r)
                    except:
                        continue
                    break

        for p in spectators['p']:
            for chem in p_match.keys():
                if p == chem:
                    try:
                        p_match.pop(p)
                    except:
                        continue
                    break
                elif p in p_match[chem]:
                    try:
                        p_match.pop(p)
                    except:
                        continue
                    break

        for p in mol_p.keys():
            if not p in p_match.keys():
                p_match[p] = []
        
        ### if molecule is not in substructure match, store in spectator dict
        for r in mol_r.keys():
            if r not in r_match.keys() and not any(mol == r for mol_l in p_match.values() for mol in mol_l):
                spectators['r'].append(r)

        for p in mol_p.keys():
            if p not in p_match.keys() and not any(mol == p for mol_l in r_match.values() for mol in mol_l):
                spectators['p'].append(p)

        ### try to make matches based on mass balance
        # create list of reaction posibilities:
        reaction_list  = generate_valid_combinations(r_match,p_match)
        if reaction_list:
            # calculate mass balance for each combination of reactants and products
            mass_balance = []
            for r,p in reaction_list:
                mass_balance.append(calculate_mass_balance(r,p))
            
            # calculate mass balance for each combination of reactants and products
            charge_balance = []
            for r,p in reaction_list:
                charge_balance.append(calculate_charge_balance(r,p))

            substructure_score = []
            for r,p in reaction_list:
                substructure_score.append(calculate_substructre_score(r,p))

            max_atom = max(x[1] for x in mass_balance)
            score = (0,0)

            for ind,r in enumerate(reaction_list):
                l_r,l_p = r
                if set(l_r) == set(l_p):
                    continue
                test_combination = score_reaction(mass_balance[ind],charge_balance[ind][1],substructure_score[ind],score,max_atom)

                if test_combination[0] == True:
                    score = (test_combination[1], mass_balance[ind][1])
                    reactants, products = r
            
            if score[1] < mass_limit:
                reactants = []
                products = []
        
        else:
            reactants = []
            products = []
        
        core_reactants_smiles = ".".join(sorted(reactants))
        core_products_smiles = ".".join(sorted(products))
        list_can_r.append(f"{core_reactants_smiles}>>{core_products_smiles}")

        for r in reactant_smiles:
            if r not in reactants and r not in spectators['r']:
                spectators['r'].append(r)
        for p in product_smiles:
            if p not in products and p not in spectators['p']:
                spectators['p'].append(p)

        # Report spectators based on first/last step or comparison adjacent steps
        if i == 0:
            dct_rogue[i] = spectators['r']
            if spectators['p']:
                rogue_molecules = spectators['p']
        else:
            if rogue_molecules == False:
                if spectators['p']:
                    rogue_molecules = spectators['p']
            else:
                r_list = []
                for r in rogue_molecules:
                    if r in spectators['r']:
                        r_list.append(r)
                if r_list:
                    dct_rogue[i] = r_list
                    
                if spectators['p']:
                    rogue_molecules = spectators['p']
                else:
                    rogue_molecules = False
        
        if i == len(list_reactions) - 1:
            dct_rogue[i+1] = spectators['p']

    dct_rogue = {k: v for k, v in dct_rogue.items() if v != []}
    return list_can_r, dct_rogue


def find_longest_matching_sequences(steps_dict):
    """Find the longest non-overlapping matching subsequences between the model and student's sequences, including shift."""
    steps = list(steps_dict.keys())
    all_match = True
    missing_steps = []
    new_dict = {}
    for x in steps:
        if len(steps_dict[x]) > 1:
            length = len(steps_dict[x])            
            if steps_dict[x][-1][0][1] == True:   # check if last step is True (for example product, or further reaction intermediate)
                new_dict[x] = ([True,True],steps_dict[x][0][1])
                miss_arr = np.arange(x+1,(x+length),1)
                missing_steps.extend(miss_arr.tolist())
                all_match = False
            elif steps_dict[x][-1][0][1] == False:
                new_dict[x] = ([True,False],steps_dict[x][0][1])
                miss_arr = np.arange(x+2,(x+length),1)
                missing_steps.extend(miss_arr.tolist())
                all_match = False
        elif steps_dict[x]:
            new_dict[x] = steps_dict[x][0]
        elif not steps_dict[x]:
            all_match = False
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
            elif new_dict[step][0] == [False,False]:
                all_match = False
    matches = list(dict.fromkeys(matches))
    matches[:] = sorted(matches, key=lambda x: (isinstance(x[0], str), x[0] if not isinstance(x[0], str) else 0))

    # Initialize tracking variables
    current_start = None
    current_end = None
    current_shift = None
    non_overlapping_matches = []

    for i, (step, shift) in enumerate(matches):
        # Helper function to check if values are consecutive integers
    
        if current_start is None:
            # Start a new sequence
            current_start = step
            current_end = step
            current_shift = shift
        # Check if current step continues the sequence (both are integers and consecutive)
        elif step != 'not_present' and current_end != 'not_present' and is_consecutive(step, current_end):
            current_end = step 
        else:
            all_match = False
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

    if not non_overlapping_matches:
        all_match = False
    
    try:
        if non_overlapping_matches[0][0] == 0 and non_overlapping_matches[0][2] != 0:
            missing_start = list(range(0,abs(non_overlapping_matches[0][2])))
            for s in missing_start:
                missing_steps.append(s)
            all_match = False
            missing_steps = sorted(list(set(missing_steps)))
    except:
        pass

    return sorted(non_overlapping_matches),missing_steps,all_match

def compare_steps(model_canonical_full, student_canonical, acid_base = None, rogue_mol = None):
    """
    Compare two lists of chemical reactions and determine how student reactions match with model reactions.
    
    Parameters:
    - model_canonical: List of model reactions in SMILES format
    - student_canonical: List of student reactions in SMILES format
    - acid_base: List of acid/base pairs that can be substituted
    
    Returns:
    - Dictionary with indices of student reactions as keys and matching information as values
    """
    if not student_canonical:
        return {0: [([False, False], 0)]}, {}, {'model_steps': len(model_canonical_full[0]), 'student_steps': 0}, {"produced":False,"last_step":False}, {}
  
    model_canonical, model_rogue = model_canonical_full
    # Get the lengths of both lists
    st_len = len(student_canonical)
    m_len = len(model_canonical)

    # Initialize the results dictionaries
    individual_comparisons = {}
    molecular_structure = {}
    product = {'produced':False,'last_step':False}
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

    reaction_found = None
    rogue_remove_r = None
    rogue_remove_p = None
    rogue_idx      = None # student match index in m_len > st_len (else)
    # Determine which approach to use based on which list is longer
    if st_len >= m_len:
        
        # Case: Student has more or equal reactions compared to model
        individual_comparisons = {key: [] for key in range(st_len)}
        molecular_structure = {key: None for key in range(st_len)}
        
        for st_idx, (st_r, st_p) in enumerate(student_reactions):
            if reaction_found:
                if isinstance(rogue_remove_r, list):
                    try:
                        for r in rogue_remove_r:
                            rogue_mol[st_idx-1].remove(r)
                    except:
                        pass
                if isinstance(rogue_remove_p, list):
                    try:
                        for p in rogue_remove_p:
                            rogue_mol[st_idx].remove(p)
                    except:
                        pass

            reaction_found = False
            reaction_candidate = False
            
            if st_r == [''] and st_p == ['']:
                individual_comparisons[st_idx] = [([False, False], 'no_chemical_transformations')]
                continue
            # Compare this student reaction with each model reaction
            for m_idx, (m_r, m_p) in enumerate(model_reactions):
                # Include reactants/products from adjacent model reactions
                try:
                    add_m_r = model_reactions[m_idx-1][1]
                except:
                    add_m_r = []
                try:
                    add_m_p = model_reactions[m_idx+1][0]
                except:
                    add_m_p = []        
                # Check reactants and products
                
                r_match_info, rogue_remove_r = check_molecules_match(st_r, m_r, add_m_r, acid_base, rogue_mol, model_rogue.get(m_idx), st_idx)
                p_match_info, rogue_remove_p = check_molecules_match(st_p, m_p, add_m_p, acid_base, rogue_mol, model_rogue.get(m_idx+1), st_idx+1)
                bool_r = r_match_info['all_match']
                bool_p = p_match_info['all_match']

                if m_len-1 == m_idx and bool_p == True:
                    if st_idx == st_len - 1:
                        product = {'produced':True,'last_step':True}
                    else:
                        product['produced'] = True
                
                # Determine the overall match status
                if bool_r and bool_p:
                    individual_comparisons[st_idx] = [([True, True], m_idx)]
                    to_be_popped = []
                    if molecular_structure[st_idx]:
                        for i_d,x in enumerate(molecular_structure[st_idx]):
                            try:
                                if x[2] == m_idx:
                                    to_be_popped.append(i_d)
                            except:
                                pass
                        for n in to_be_popped:
                            molecular_structure[st_idx].pop(n)
                    reaction_found = True     
                    break
                
                elif bool_r:
                    individual_comparisons[st_idx].append(([True, False], m_idx))
                    if False in p_match_info['individual_matches'] and p_match_info['unique_molecules'] and p_match_info['unique_from_model']:
                        differences,mol_rogue = analyze_molecular_differences(
                            p_match_info['unique_molecules'], 
                            p_match_info['unique_from_model'],
                            rogue_mol,
                            st_idx+1,)

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
                   
                    if False in r_match_info['individual_matches'] and r_match_info['unique_molecules'] and r_match_info['unique_from_model']:
                        
                        if molecular_structure[st_idx] != None:
                            pass
                        else:
                            differences,mol_rogue = analyze_molecular_differences(
                                r_match_info['unique_molecules'], 
                                r_match_info['unique_from_model'],
                                rogue_mol,
                                st_idx)

                            if all(x == None for x in differences):
                                pass    
                            else:
                                molecular_structure[st_idx] = differences
                    continue
                
                elif True in r_match_info['individual_matches']:
                    reaction_candidate = True
                    candidate_index = st_idx
                    r_match_cand = r_match_info
                    m_cand_idx = m_idx

            # If no match found for this student reaction
            if not reaction_found and reaction_candidate == True:
                individual_comparisons[candidate_index].append(([False, False], m_cand_idx))
                if molecular_structure[candidate_index] != None:
                    pass
                else:
                    differences,mol_rogue = analyze_molecular_differences(
                        r_match_cand['unique_molecules'], 
                        r_match_cand['unique_from_model'],
                        rogue_mol,
                        candidate_index)

                    if all(x == None for x in differences):
                        pass    
                    else:
                        molecular_structure[candidate_index] = differences

            if not reaction_found and not reaction_candidate:
                individual_comparisons[st_idx].append(([False, False], 'not_present'))
        if reaction_found:
            if isinstance(rogue_remove_r, list):
                for r in rogue_remove_r:
                    try:
                        rogue_mol[st_idx].remove(r)
                    except:
                        pass
            if isinstance(rogue_remove_p, list):
                for p in rogue_remove_p:
                    try:
                        rogue_mol[st_idx+1].remove(p)
                    except:
                        pass
    else:
        # Case: Model has more reactions than student
        individual_comparisons = {key: [] for key in range(st_len+1)}
        molecular_structure = {key: None for key in range(st_len+1)}
        # Consider whether the final compound ensamble matches the model
        final_step = []
        try:
            final_step.append(student_reactions[0][-1][0])
            final_step.append(rogue_mol.get(st_len)[0])
        except:
            pass

        student_reactions.append((final_step,[]))
        rogue_mol.pop(st_len)

        for m_idx, (m_r, m_p) in enumerate(model_reactions):
            if reaction_found:
                if isinstance(rogue_remove_r, list):
                    try:
                        for r in rogue_remove_r:
                            rogue_mol[st_idx].remove(r)
                    except:
                        pass
                if isinstance(rogue_remove_p, list):
                    try:
                        for p in rogue_remove_p:
                            rogue_mol[st_idx+1].remove(p)
                    except: 
                        pass
            
            reaction_found = False
            reaction_candidate = False
            for st_idx, (st_r, st_p) in enumerate(student_reactions):
                # Skip if we already found a perfect match for this student reaction
                if st_idx in individual_comparisons and individual_comparisons[st_idx] and individual_comparisons[st_idx][0][0] == [True, True]:
                    continue
                # Include reactants/products from adjacent model reactions
                try:
                    add_m_r = model_reactions[m_idx-1][1]
                except:
                    add_m_r = []
                try:
                    add_m_p = model_reactions[m_idx+1][0]
                except:
                    add_m_p = []        
                # Check reactants and products
                #rogue_student
                r_match_info, rogue_remove_r = check_molecules_match(st_r, m_r, add_m_r, acid_base, rogue_mol, model_rogue.get(m_idx), st_idx)
                p_match_info, rogue_remove_p = check_molecules_match(st_p, m_p, add_m_p, acid_base, rogue_mol, model_rogue.get(m_idx+1), st_idx+1)
                bool_r = r_match_info['all_match']
                bool_p = p_match_info['all_match']
                if m_len-1 == m_idx and bool_p == True:
                    if st_idx == st_len - 1:
                        product = {'produced':True,'last_step':True}
                    else:
                        product['produced'] = True
                    
                # Determine the overall match status
                if bool_r and bool_p:
                    individual_comparisons[st_idx] = [([True, True], m_idx)]
                    to_be_popped = []
                    if molecular_structure[st_idx]:
                        for i_d,x in enumerate(molecular_structure[st_idx]):
                            try:
                                if x[2] == m_idx:
                                    to_be_popped.append(i_d)
                            except:
                                pass
                        for n in to_be_popped:
                            molecular_structure[st_idx].pop(n)
                    reaction_found = True
                    rogue_index = st_idx
                    break
                elif bool_r:
                    individual_comparisons[st_idx].append(([True, False], m_idx))
                    if True in p_match_info['individual_matches'] and False in p_match_info['individual_matches'] and p_match_info['unique_molecules'] and p_match_info['unique_from_model']:
                        differences, rogue_mol = analyze_molecular_differences(
                            p_match_info['unique_molecules'], 
                            p_match_info['unique_from_model'],
                            rogue_mol,
                            st_idx+1)

                        if all(x == None for x in differences):
                            pass    
                        else:
                            molecular_structure[st_idx+1] = differences

                    rogue_index = st_idx
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
                            differences,mol_rogue = analyze_molecular_differences(
                                r_match_info['unique_molecules'], 
                                r_match_info['unique_from_model'],
                                rogue_mol,
                                st_idx)
                            if all(x == None for x in differences):
                                pass    
                            else:
                                molecular_structure[st_idx] = differences 
                    rogue_index = st_idx
                    continue

                elif not individual_comparisons[st_idx] and True in r_match_info['individual_matches']:
                    reaction_candidate = True
                    candidate_index = st_idx
                    r_match_cand = r_match_info
                    m_cand_idx = m_idx
                    rogue_index = st_idx
            # If no match found for this student reaction
            if not reaction_found and reaction_candidate == True:
                individual_comparisons[candidate_index].append(([False, False], m_cand_idx))
                if molecular_structure[candidate_index] != None:
                    pass
                else:
                    differences,mol_rogue = analyze_molecular_differences(
                        r_match_cand['unique_molecules'], 
                        r_match_cand['unique_from_model'],
                        rogue_mol,
                        candidate_index)
                    if all(x == None for x in differences):
                        pass    
                    else:
                        molecular_structure[candidate_index] = differences
                    rogue_index = st_idx

            if not reaction_found and not reaction_candidate:
                if all_lists_empty(individual_comparisons):
                    individual_comparisons[0].append(([False, False], "not_present"))
                else:
                    ind_add = max([k for k, v in individual_comparisons.items() if v]) # find key in dict with a not empty list
                    individual_comparisons[ind_add].append(([False, False], "not_present"))
        if reaction_found:
            if isinstance(rogue_remove_r, list):
                for r in rogue_remove_r:
                    try:
                        rogue_mol[rogue_index].remove(r)
                    except:
                        pass
            if isinstance(rogue_remove_p, list):
                for p in rogue_remove_p:
                    try:
                        rogue_mol[rogue_index+1].remove(p)
                    except:
                        pass
                
    # Clean up the molecular structure dictionary

    for r in individual_comparisons.keys():
        if individual_comparisons[r] and 'no_chemical_transformations' in individual_comparisons[r][0]:
            molecular_structure.pop(r)
    molecular_structure = {k: v for k, v in molecular_structure.items() if v is not None}
    
    # Add general information
    general = {
        "model_steps": m_len,
        "student_steps": st_len,
    }
    ### remove empty lists as values from individual_comparisons
    individual_comparisons = {k: v for k, v in individual_comparisons.items() if v}
    rogue_mol = {k: v for k, v in rogue_mol.items() if v != []}
    
    return individual_comparisons, molecular_structure, general, product, rogue_mol


def check_molecules_match(molecules1, molecules2, add_molecules2, acid_base, rogue_student, rogue_model, index):
    """
    Check how molecules from list1 match with molecules from list2, considering acid/base substitutions.
    
    Returns a dictionary with match information.
    """
    individual_matches = []
    matching_molecules = []
    unique_molecules = []
    missing_molecules = []

    # Check matches from molecules1 to molecules2
    for mol1 in molecules1:
        if mol1 in molecules2 or mol1 in add_molecules2:
            individual_matches.append(True)
            matching_molecules.append(mol1)
        elif any(mol1 in acid for acid in acid_base[0]) and any(acid in molecules2 for acid in acid_base[0]) or any(mol1 in acid for acid in acid_base[0]) and any(acid in add_molecules2 for acid in acid_base[0]):
            individual_matches.append(True)
            matching_molecules.append(mol1)
        elif any(mol1 in base for base in acid_base[1]) and any(base in molecules2 for base in acid_base[1]) or any(mol1 in base for base in acid_base[1]) and any(base in add_molecules2 for base in acid_base[1]):
            individual_matches.append(True)
            matching_molecules.append(mol1)
        else:
            individual_matches.append(False)
            unique_molecules.append(mol1)
    
    to_remove = []

    # Check for molecules in molecules2 that don't have matches in molecules1
    for mol2 in molecules2:
        found_match = False
        # Direct match
        try:
            matches_acid = set(rogue_student.get(index)) & set(acid_base[0])
        except:
            matches_acid = []
        try:
            matches_base = set(rogue_student.get(index)) & set(acid_base[1])
        except:
            matches_base = []
        if mol2 in molecules1:
            found_match = True
        elif rogue_student.get(index) and mol2 in rogue_student.get(index):
            individual_matches.append(True)
            matching_molecules.append(mol2)
            to_remove.append(mol2)
            found_match = True
        # Check rogue_student in acid and mol2 in acid.
        elif matches_acid and any(mol2 in acid for acid in acid_base[0]):
            r = elem = next(iter(matches_base), None)
            individual_matches.append(True)
            matching_molecules.append(r)
            to_remove.append(r)
            found_match = True
        elif matches_base and any(mol2 in base for base in acid_base[1]):
            r = elem = next(iter(matches_base), None)
            individual_matches.append(True)
            matching_molecules.append(r)
            to_remove.append(r)
            found_match = True
        # Check acid substitutions
        elif any(mol2 in acid for acid in acid_base[0]):
            to_remove.append(mol2)
            found_match = True
        # Check base substitutions  
        elif any(mol2 in base for base in acid_base[1]):
            found_match = True

        if not found_match:
            missing_molecules.append(mol2)

    if rogue_student.get(index) and rogue_model:
        try:
            for s in rogue_student.get(index):
                remove = False
                for m in rogue_model:
                    remove = mols_equal(Chem.MolFromSmiles(s),Chem.MolFromSmiles(m))
                    if remove == True:
                        to_remove.append(s)
                        break
        except:
            pass

    unique_from_model = set(molecules2) - set(matching_molecules)
    
    # all_match is True only if all molecules1 match AND no molecules are missing from molecules2
    all_match = all(individual_matches) and len(missing_molecules) == 0

    return {
        'individual_matches':   individual_matches,
        'all_match'         :   all_match,
        'matching_molecules':   matching_molecules,
        'unique_molecules'  :   unique_molecules,
        'missing_molecules' :   missing_molecules,
        'unique_from_model' :   unique_from_model
    }, to_remove


def analyze_molecular_differences(unique_molecules, unique_from_model, rogue, index):
    """
    Analyze differences between molecules that don't match directly.
    """
    if not unique_from_model:
        # If there's no matching molecule in the model, report it
        return ['molecule_not_present_in_model', sum(1 for atom in Chem.MolFromSmiles(unique_molecules[0]).GetAtoms() if atom.GetSymbol() == 'C')]
    
    to_remove = []
    try:
        for r in rogue[index]:
            r_mol = Chem.MolFromSmiles(r)
            for m in unique_from_model:
                m_mol = Chem.MolFromSmiles(m)
                if r_mol.GetNumAtoms() == m_mol.GetNumAtoms() and r_mol.HasSubstructMatch(m_mol) and m_mol.HasSubstructMatch(r_mol):
                    to_remove.append(r)
        for r in to_remove:
            rogue = rogue[index].remove(r)
    except:
        pass

    
    differences = []
    for mol in unique_molecules:
        comparison = []
        for model_mol in unique_from_model:
            comparison.append(compare_molecules(mol, model_mol))
        # Find the best match based on priority
 
        pri = 8
        diff = None
        for z in comparison:
            if z[0].endswith('_1') and pri > 0:
                diff = [z[0].strip('_1'), z[1],10]
                pri = 0
            elif z[0].endswith('_2') and pri > 1:
                diff = [z[0].strip('_2'), z[1],9]
                pri = 1
            elif z[0].endswith('_3') and pri > 2:
                diff = [z[0].strip('_3'), z[1],8]
                pri = 2
            elif z[0].endswith('_4') and pri > 3:
                diff = [z[0].strip('_4'), z[1],7]
                pri = 3
            elif z[0].endswith('_5') and pri > 4:
                diff = [z[0].strip('_5'), z[1],6]
                pri = 4
            elif z[0].endswith('_6') and pri > 5:
                diff = [z[0].strip('_6'), z[1],5]
                pri = 5    
            elif z[0].endswith('_7') and pri > 6:
                diff = [z[0].strip('_7'), z[1],4]
                pri = 6
            elif z[0].endswith('_0') and pri > 7:
                diff = [z[0].strip('_0'), z[1],3]
                pri = 7
        
        differences.append(diff)

    return differences,rogue

def compare_reactions(model_list, student_list, acid_base):
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
    student_canonical = canonicalize_reaction_list(student_list)

    index_comp = False
    rogue_mols = student_canonical[1]
    if isinstance(model_list, list) and all(isinstance(item, list) for item in model_list):
        lst_comparisons = []
        lst_m_canonical = []
        for m_lst in model_list:
            model_canonical = canonicalize_reaction_list(m_lst)
            lst_m_canonical.append(model_canonical[0])
            rogue_mols_copy = copy.deepcopy(rogue_mols)
            ind_comp = compare_steps(model_canonical,student_canonical[0],acid_base,rogue_mols_copy)
            lst_comparisons.append(ind_comp)
                # Count keys in the first dictionary of each tuple
        key_counts = [len(comparison[0]) for comparison in lst_comparisons]
        
        # Check if there's a tie for maximum number of keys
        max_count = max(key_counts)
        has_tie = key_counts.count(max_count) > 1 and key_counts[0] != 0
        bool_t_f = []
        weight_struc = []

        for x in lst_comparisons:
            true_count = 0
            false_count = 0
            w = 0
            for values in x[0].values():
                for condition_list, _ in values:
                    true_count += condition_list.count(True)
                    false_count += condition_list.count(False)
            bool_t_f.append(true_count - false_count)
            for struc in x[1].keys():
                for weight in x[1][struc]:
                    w += weight[2]
            if w == 0:
                w = 1000
            weight_struc.append(w)

        if has_tie and bool_t_f[0] != bool_t_f[1]:
            index_comp, val = max(enumerate(bool_t_f), key=lambda x: x[1])

            # Select the one with max True
            individual_comparisons = lst_comparisons[index_comp]
        elif has_tie and weight_struc[0] != weight_struc[1]:
            index_comp,val = max(enumerate(weight_struc), key=lambda x: x[1])
            individual_comparisons = lst_comparisons[index_comp]
            
        elif has_tie and bool_t_f[0] == bool_t_f[1]:
            index_comp = key_counts.index(max_count)
            individual_comparisons = lst_comparisons[index_comp]
        else:
            index_comp, individual_comparisons = max(enumerate(lst_comparisons), key=lambda x: count_trues(x[1]))
        
        m_canonical = (lst_m_canonical[index_comp], {})
        

    else:   
        m_canonical = canonicalize_reaction_list(model_list)
        individual_comparisons = compare_steps(m_canonical, student_canonical[0],acid_base,student_canonical[1])

    matching_subsequences = find_longest_matching_sequences(individual_comparisons[0])
    
    if not matching_subsequences:
        return {"individual_steps": individual_comparisons[0],"matching_sequences": [],
                "missing_steps":[],"molecular_structures":individual_comparisons[1],"index_resonance": index_comp,"all_steps_match":False, "product": individual_comparisons[3],'reactions_model':m_canonical,'reactions_student':student_canonical[0],'rogue_species':individual_comparisons[4]}
    else:
        return {"individual_steps": individual_comparisons[0],"matching_sequences": matching_subsequences[0],
                "missing_steps":matching_subsequences[1],"molecular_structures":individual_comparisons[1],"index_resonance": index_comp,"all_steps_match":matching_subsequences[2], "product": individual_comparisons[3],'reactions_model':m_canonical,'reactions_student':student_canonical[0],'rogue_species':individual_comparisons[4]}