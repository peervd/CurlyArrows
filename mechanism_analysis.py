import numpy as np
from rdkit import Chem
from mechanism_comparison import compare_reactions, canonicalize_reaction
from exersices_ORC import load_acid_base
from collections import Counter
from difflib import SequenceMatcher

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def unique_strings(list_of_lists):
    # Flatten the list and count occurrences
    all_strings = [item for sublist in list_of_lists for item in set(sublist)]
    counts = Counter(all_strings)
    
    # Keep only elements that appear in exactly one sublist
    unique_list = [item for sublist in list_of_lists for item in sublist if counts[item] == 1]
    
    return unique_list

def common_substring(strings):
    if len(strings[0]) > 6 and len(strings[1]) > 6:
        matcher = SequenceMatcher(None, strings[0], strings[1])
        for match in matcher.get_matching_blocks():
            if match.size >= 5:
                return True
    return False

def flatten(xss):
    return [x for xs in xss for x in xs]

def custom_sort(item):
    if item == 'reactants':
        return (0, 0)
    elif item == 'products':
        return (1, 0)
    elif item.startswith('intermediates_'):
        parts = item.split('_')
        number = int(parts[1]) 
        return (2, number)
    return (3, 0)  

def filter_variations(lst):
    base_items = [item for item in lst if not item.endswith(("_a", "_b"))]
    grouped_pairs = {}

    # Group _a and _b pairs by number
    for item in lst:
        if item.endswith("_a") or item.endswith("_b"):
            key = item.rsplit("_", 1)[0]  # Extract base identifier (e.g., "intermediates_1")
            grouped_pairs.setdefault(key, []).append(item)

    # Ensure each key has exactly two elements (_a and _b)
    grouped_pairs = {k: v for k, v in grouped_pairs.items() if len(v) == 2}
    variations = []

    # Generate variations by keeping one element from each pair
    keys = list(grouped_pairs.keys())
    for i in range(2 ** len(keys)):  # Iterate over all possible 2^n combinations
        variation = base_items[:]
        for j, key in enumerate(keys):
            variation.append(grouped_pairs[key][(i >> j) & 1])  # Choose _a or _b based on binary representation
        
        sorted_data = sorted(variation, key=custom_sort)
        variations.append(sorted_data)
    return variations

def generate_reactions(reaction_dict, type = 'model'):
    reaction_list = []
    steps = list(reaction_dict.keys())
    n_steps = len(steps)
    resonance = any(step.endswith('a') for step in steps)
    variations_steps = False
    # Extract molecules for each step
    mols = {s: list(reaction_dict[s]['molecules'].values()) for s in steps}
    if n_steps > 2:
        if resonance == True and type == 'model':
            variations_steps = filter_variations(steps)
            for i, v in enumerate(variations_steps, 1):
                sub_reaction = []
                n_steps_v = len(variations_steps[i-1])
                for it, step in enumerate(v[:-1]):
                    
                    if it == 0:
                        reactants = mols[v[0]]
                        products = mols[v[2]]
                    elif it == (n_steps_v - 2):
                        reactants = mols[v[it + 1]]
                        products = mols[v[1]]
                    else:
                        reactants = mols[v[it + 1]]
                        products = mols[v[it + 2]]
        
                    # Identify unique molecules in reactants and products
                    reactant_set = set(reactants)
                    product_set = set(products)
                    
                    reactants_str = ".".join(m for m in reactant_set if m not in product_set)
                    products_str = ".".join(m for m in product_set if m not in reactant_set)
        
                    reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                    reaction = canonicalize_reaction(reaction)
                    sub_reaction.append(reaction)
                reaction_list.append(sub_reaction)

        elif resonance == True and type == 'student':
            filtered_steps = [step for step in steps if not step.endswith('b')]
            n_steps_v = len(filtered_steps)
            for i, step in enumerate(filtered_steps[:-1]):
                if i == 0:
                    reactants = mols[filtered_steps[0]]
                    products = mols[filtered_steps[2]]
                elif i == (n_steps_v - 2):
                    reactants = mols[filtered_steps[i + 1]]
                    products = mols[filtered_steps[1]]
                else:
                    reactants = mols[filtered_steps[i + 1]]
                    products = mols[filtered_steps[i + 2]]
    
                # Identify unique molecules in reactants and products
                reactant_set = set(reactants)
                product_set = set(products)
                
                reactants_str = ".".join(m for m in reactant_set if m not in product_set)
                products_str = ".".join(m for m in product_set if m not in reactant_set)
    
                reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                reaction = canonicalize_reaction(reaction)
                reaction_list.append(reaction)
            
        elif resonance == False:     
            for i, step in enumerate(steps[:-1]):
                if i == 0:
                    reactants = mols[steps[0]]
                    products = mols[steps[2]]
                elif i == (n_steps - 2):
                    reactants = mols[steps[i + 1]]
                    products = mols[steps[1]]
                else:
                    reactants = mols[steps[i + 1]]
                    products = mols[steps[i + 2]]
    
                # Identify unique molecules in reactants and products
                reactant_set = set(reactants)
                product_set = set(products)
                
                reactants_str = ".".join(m for m in reactant_set if m not in product_set)
                products_str = ".".join(m for m in product_set if m not in reactant_set)
    
                reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                reaction = canonicalize_reaction(reaction)
                reaction_list.append(reaction)
                
    elif n_steps == 2:
        reactants = mols[steps[0]]
        products = mols[steps[1]]

        # Identify unique molecules in reactants and products
        reactant_set = set(reactants)
        product_set = set(products)
        
        reactants_str = ".".join(m for m in reactant_set if m not in product_set)
        products_str = ".".join(m for m in product_set if m not in reactant_set)

        reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
        reaction = canonicalize_reaction(reaction)
        reaction_list.append(reaction)        
    elif n_steps < 2:
        reaction = "No reaction takes place"
    if variations_steps == False:
        return reaction_list,None
    else:
        return reaction_list,variations_steps


def reaction_rules(arrows,reaction):
    ''' Assesses whether a reaction step can be identified by a common
    reaction transformation: proton transfer, nucleophilic attck, loss of 
    leaving group, rearrangement, redox'''

    reaction = reaction.split('>>')
    reactants = reaction[0]
    products =  reaction[1]
    r = reactants.split('.')
    p = products.split('.')
    r_chem = []
    for x in r:
        mol = Chem.MolFromSmiles(x,sanitize=False)
        r_chem.append(mol)
    p_chem = []
    for x in p:
        mol = Chem.MolFromSmiles(x,sanitize=False)
        p_chem.append(mol)
    
    reaction_type = []

    for a in arrows:

        ''' discriminate proton transfer from rearrangemnt '''
        if 'H' in a[0] or 'H' in a[1]:
            reaction_type.append('proton transfer')
            continue

        ''' identify nucleophilic attack '''
        if 'H' not in a[0] and 'H' not in a[1]:
            if a[4] == a[5] and common_substring(a[2:4]) == False:
                reaction_type.append('nucleophilic attack')
                continue
            elif a[4] != a[5]:
                reaction_type.append('nucleophilic attack')
                continue 
            
        ''' identify loss of leaving group / rearrangemnt'''
        if 'H' not in a[0] and 'H' not in a[1] and a[2] == a[3]:
            if max(s.GetNumAtoms() for s in r_chem) > max(s.GetNumAtoms() for s in p_chem):
                reaction_type.append('loss of leaving group')
            elif max(s.GetNumAtoms() for s in r_chem) == max(s.GetNumAtoms() for s in p_chem):
                reaction_type.append('rearrangement')
            else:
                reaction_type.append('internal electron movement')  

    return reaction_type


def compare_arrows(model_arrows, student_arrows,model_reaction,student_reaction):
    if not model_arrows or not student_arrows:  # Handle empty input
        if model_arrows:
            model_arr = []
            for x in range(0,len(model_arrows)):
                arrows = [(model_arrows[x].get('start_arrow')),
                (model_arrows[x].get('end_arrow')),
                (model_arrows[x].get('substructure_start')),
                (model_arrows[x].get('substructure_end')),
                (model_arrows[x].get('start_molecule')),
                (model_arrows[x].get('end_molecule'))]
                model_arr.append(arrows)
            return [reaction_rules(model_arr,model_reaction),[False,False]]   
        else:
            return [[False, False],[False,False]]

    largest_index, largest_value = max(enumerate([len(model_arrows),len(student_arrows)]), key=lambda x: x[1])

    if largest_index == 0:
        long_mech = model_arrows
        short_mech = student_arrows
        student = 'short_mechanism'
    elif largest_index == 1:
        long_mech = student_arrows
        short_mech = model_arrows
        student = 'long_mechanism'

    long_arr = []
    for x in range(0,len(long_mech)):
        arrows = [(long_mech[x].get('start_arrow')),
        (long_mech[x].get('end_arrow')),
        (long_mech[x].get('substructure_start')),
        (long_mech[x].get('substructure_end')),
        (long_mech[x].get('start_molecule')),
        (long_mech[x].get('end_molecule'))]
        long_arr.append(arrows)

    short_arr = []
    for x in range(0,len(short_mech)):
        arrows = [(short_mech[x].get('start_arrow')),
        (short_mech[x].get('end_arrow')),
        (short_mech[x].get('substructure_start')),
        (short_mech[x].get('substructure_end')),
        (short_mech[x].get('start_molecule')),
        (short_mech[x].get('end_molecule'))]
        short_arr.append(arrows)

    matches = []
    max_it = len(short_arr) - 1
    for x in long_arr:
        app = False
        for i,y in enumerate(short_arr):
            if x[0:4] == y[0:4]:
                matches.append([ True,True ])
                app = True
                break
            elif (x[0],x[2]) == (y[0],y[2]):
                matches.append([ True,False ])
                app = True
                break
            elif [x[1],x[3]] == [y[1],y[3]]:
                matches.append([ False,True ])
                app = True
                break
            elif x == [y[1],y[0],y[3],y[2]]:
                matches.append([ 'arrow_flipped',False ])
                app = True
                break
            elif app == False and max_it == i:
                matches.append([ False,False ])
 
    if any(x == False for y in matches for x in y) and not any(x == 'arrow_flipped' for y in matches for x in y):
        if student == 'short_mechanism':
            model_operation = reaction_rules(long_arr,model_reaction)
            student_operation = reaction_rules(short_arr,student_reaction)
            matches = [model_operation,student_operation]
            if matches[0] == matches[1]:
                matches = [[True,True],[True,True]]

        elif student == 'long_mechanism':
            model_operation = reaction_rules(short_arr,model_reaction)
            student_operation = reaction_rules(long_arr,student_reaction)
            matches = [model_operation,student_operation]
            if matches[0] == matches[1]:
                matches = [[True,True],[True,True]]
    return matches

def check_resonance(model, student, steps):
    # If no keys in 'student' end with '_a', return 'no_resonance'
    if not any(key.endswith('_a') for key in student):
        return ['no_resonance']

    resonance = []
    
    for key in student:
        if key.endswith('_b'):
            i = key.split('_')
            step_index = float(i[1])
            for x in steps['matching_sequences']:
                if (x[0] - x[2]) <= step_index <= (x[1] - x[2]):
                    m_int = f"intermediates_{int(step_index) - x[2]}"

                    # Find all relevant model keys
                    m_res = [k for k in model if m_int in k]
                    
                    # Extract molecules from model steps
                    step_molecules = [list(model[step]['molecules'].values()) for step in m_res]
                    
                    # Find unique resonance molecules
                    resonance_molecules = unique_strings(step_molecules)

            # Get student's resonance molecules
            resonance_student = list(student[key]['molecules'].values())
            # Check if student's resonance molecule matches any valid ones

            if resonance_molecules and resonance_student[0] in resonance_molecules:
                resonance.append([True,key])
            else:
                resonance.append([False,key])

    return resonance

def individual_steps(model,student,resonance):
    
    ''' 
    The student answer of the reactin mechanism is assessed on a global level.
    First, individual reaction steps from both the model answer and student
    answer are generated as reaction transformations in SMILES. Next, the 
    individual steps between the model and student answer are compared. If 
    all steps between the student and model match, the mechanism is correct
    on a global level. Otherwise the inconsistencies are returned in a 
    dictionary.
    
    '''
    model_reactions = generate_reactions(model, type = 'model')
    student_reactions = generate_reactions(student, type = 'student')
    
    steps = compare_reactions(model_reactions[0],student_reactions[0])
    
    s_keys = steps['individual_steps'].keys()
    
    steps_bool = []
    for x in s_keys:
        steps_bool.append(steps['individual_steps'][x][0])
    
    exersice = True
    if all(all(not item for item in sublist) for sublist in steps_bool): 
        
        model_m = []
        for x in list(model['reactants']['molecules'].keys()):
            model_m.append(model['reactants']['molecules'][x].replace('[H]',''))
        student_m = []
        for x in list(student['reactants']['molecules'].keys()):
            student_m.append(student['reactants']['molecules'][x].replace('[H]',''))
        correct_exercise = []
        
        for x in model_m:
            if x in student_m and x not in flatten(load_acid_base('list_acid_base.csv')):
                correct_exercise.append(True)
            else:
                correct_exercise.append(False)

        if all(x == False for x in correct_exercise):
            exersice = False
    ### check whether resonance structures are present
    resonance = []
    presence_resonance = []
    if exersice == True:
        presence_resonance = check_resonance(model,student,steps)
        
        if resonance == 'yes' or resonance == 'Yes' or resonance == 'YES':
            resonance = True
        else: 
            resonance = False

        
    if steps['index_resonance'] is False:
        return steps,model_reactions[0],student_reactions[0],exersice,list(model.keys()),{'resonance':presence_resonance,'resonance_present':resonance}
    else:
        return steps,model_reactions[0][steps['index_resonance']],student_reactions[0],exersice,model_reactions[1][steps['index_resonance']],{'resonance':presence_resonance,'resonance_present':resonance}

def reaction_transformations(model,student,steps):
    seq = steps[0]['matching_sequences']
    transformations = []
    
    first_step = False
    model_keys = steps[4]
    student_keys = [item for item in list(student.keys()) if not item.endswith(("_b"))]
    len_steps = len(student_keys) - 1
    
    if seq[0][0] != 0:
        m = model['reactants'].get('arrows')
        s = student['reactants'].get('arrows')
        check = compare_arrows(m,s,steps[1][0],steps[2][0])
        transformations.append(check)
        first_step = True
    
    for series in seq:
        for step in np.arange(series[0],series[1]+2, 1):
            if step == 0 and first_step == False:
                m = model['reactants']['arrows']
                if series[2] == 0:
                    s = student['reactants']['arrows']
                else:
                    comp = student_keys[step+1+series[2]]
                    s = student[comp]['arrows']
                check = compare_arrows(m,s,steps[1][step-series[2]],steps[2][step])
                transformations.append(check)
                if len(list(model.keys())) < 3:
                    break
            
            elif step > 0 and step < len_steps: ### if error: series[2] changed to series [1]
 
                comp_m = model_keys[step-(series[2])+1]
                if series[2] == 0:
                    comp_s = student_keys[step+1]
                elif series[2] != 0:
                    comp_s = student_keys[step+1]

                m = model[comp_m]['arrows']
                s = student[comp_s]['arrows']

                check = compare_arrows(m,s,steps[1][step-series[2]],steps[2][step])
                transformations.append(check)
            elif step >= len_steps:
                break
            
  
    return transformations
            
