import numpy as np
from rdkit import Chem
from mechanism_comparison import compare_reactions, canonicalize_reaction
from exersices_ORC import load_acid_base
from helper_functions import unique_strings, common_substring, flatten, filter_variations_main_res, filter_variations_res, custom_reaction_sort, _split_reaction

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def generate_reactions(reaction_dict, reaction_type = 'model'):
    reaction_list = []
    steps = sorted(list(reaction_dict.keys()),key = custom_reaction_sort)
    n_steps = len(steps)
    resonance = any(step.endswith('ra') for step in steps)
    main_var = any(step.endswith('ma') for step in steps)
    variations_steps = False
    # Extract molecules for each step
    mols = {s: list(reaction_dict[s]['molecules'].values()) for s in steps}
    if n_steps > 2:
        if resonance == True and main_var == True and reaction_type == 'model':
            variations_steps = filter_variations_main_res(steps)
            for i, v in enumerate(variations_steps, 1):
                sub_reaction = []
                n_steps_v = len(variations_steps[i-1])
                for it, step in enumerate(v[:-1]):
                    reactants = mols[v[it]]
                    products = mols[v[it+1]]
                    # Identify unique molecules in reactants and products
                    reactant_set = set(reactants)
                    product_set = set(products)
                    
                    reactants_str = ".".join(m for m in reactant_set if m not in product_set)
                    products_str = ".".join(m for m in product_set if m not in reactant_set)
        
                    reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                    reaction = canonicalize_reaction(reaction)
                    sub_reaction.append(reaction)
                reaction_list.append(sub_reaction)

                
        elif resonance == True and reaction_type == 'model':
            variations_steps = filter_variations_res(steps)
            for i, v in enumerate(variations_steps, 1):
                sub_reaction = []
                n_steps_v = len(variations_steps[i-1])
                for it, step in enumerate(v[:-1]):
                    reactants = mols[v[it]]
                    products = mols[v[it+1]]
        
                    # Identify unique molecules in reactants and products
                    reactant_set = set(reactants)
                    product_set = set(products)
                    
                    reactants_str = ".".join(m for m in reactant_set if m not in product_set)
                    products_str = ".".join(m for m in product_set if m not in reactant_set)
        
                    reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                    reaction = canonicalize_reaction(reaction)
                    sub_reaction.append(reaction)
                reaction_list.append(sub_reaction)

        elif resonance == True and reaction_type == 'student':
            filtered_steps = [step for step in steps if not step.endswith('rb')]
            n_steps_v = len(filtered_steps)
            for i, step in enumerate(filtered_steps[:-1]):

                reactants = mols[filtered_steps[i]]
                products = mols[filtered_steps[i+1]]
              
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

                reactants = mols[steps[i]]
                products = mols[steps[i+1]]

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
    r,p = _split_reaction(reaction)
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
        ''' discriminate proton transfer from nuc attack '''

        if 'H' in a[0] or 'H' in a[1]:
            if 'H-' not in a[2]:
                reaction_type.append('proton transfer')
                continue
            elif 'H-' == a[2]:
                reaction_type.append('nucleophilic attack')
                continue
            
            ''' identify nucleophilic attack '''
        else:
            if a[4] == a[5]:
                if max(s.GetNumAtoms() for s in r_chem) > max(s.GetNumAtoms() for s in p_chem):
                    reaction_type.append('loss of leaving group')
                    continue
                elif max(s.GetNumAtoms() for s in r_chem) == max(s.GetNumAtoms() for s in p_chem):
                    reaction_type.append('rearrangement')
                    continue

            elif a[4] != a[5]:
                reaction_type.append('nucleophilic attack')
                continue 

    return reaction_type

def reaction_arrow_check(steps,arrows):
    step_numbers = steps.keys()
    wrong_arrows = {}
    for s in step_numbers:
        if steps[s][0][0] == True:
            m_step = steps[s][1]
            if m_step == 'not_present':
                continue
            if arrows['s_arrows'][s] != arrows['m_arrows'][m_step] and arrows['m_arrows'][m_step] == 'synthetic':
                wrong_arrows[int(s)+1] = False
    return wrong_arrows

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
                (model_arrows[x].get('end_molecule')),
                (model_arrows[x].get('electrons'))]
                model_arr.append(arrows)
            return ([reaction_rules(model_arr,model_reaction),[False,False]],False,False)   
        else:
            return ([[False, False],[False,False]],False,False)

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
        (long_mech[x].get('end_molecule')),
        (long_mech[x].get('electrons'))]
        long_arr.append(arrows)

    short_arr = []
    for x in range(0,len(short_mech)):
        arrows = [(short_mech[x].get('start_arrow')),
        (short_mech[x].get('end_arrow')),
        (short_mech[x].get('substructure_start')),
        (short_mech[x].get('substructure_end')),
        (short_mech[x].get('start_molecule')),
        (short_mech[x].get('end_molecule')),
        (short_mech[x].get('electrons'))]
        short_arr.append(arrows)
    
    incorrect_arrow = False
    step_correctness = False
    matches = []
    max_it = len(short_arr) - 1
    
    for x in long_arr:
        app = False
        for i,y in enumerate(short_arr):
            if [x[0],x[1],x[2],x[3]] == [y[0],y[1],y[2],y[3]]:
                step_correctness = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                matches.append([ True,True ])
                app = True
                break
            elif (x[0],x[2]) == (y[0],y[2]):
                if x[6] != y[6]:
                    incorrect_arrow = True
                matches.append([ True,False ])
                app = True
                break
            elif [x[1],x[3]] == [y[1],y[3]]:
                if x[6] != y[6]:
                    incorrect_arrow = True
                matches.append([ False,True ])
                app = True
                break
            elif [x[0],x[1],x[2],x[3]] == [y[1],y[0],y[3],y[2]]:
                if x[6] != y[6]:
                    incorrect_arrow = True
                matches.append([ 'arrow_flipped',False ])
                app = True
                break
            elif app == False and max_it == i:
                matches.append([ False,False ])
                
                
    if student_reaction and any(x == False for y in matches for x in y) and not any(x == 'arrow_flipped' for y in matches for x in y):
        if student == 'short_mechanism':
            model_operation = reaction_rules(long_arr,model_reaction)
            student_operation = reaction_rules(short_arr,student_reaction)
            matches = [model_operation,student_operation]
            if matches[0] == matches[1] or all(item == 'proton transfer' for y in matches for item in y):
                matches = [[True,True],[True,True]]
                step_correctness = True

        elif student == 'long_mechanism':
            model_operation = reaction_rules(short_arr,model_reaction)
            student_operation = reaction_rules(long_arr,student_reaction)
            matches = [model_operation,student_operation]
            if matches[0] == matches[1] or all(item == 'proton transfer' for y in matches for item in y):
                matches = [[True,True],[True,True]]
                step_correctness = True
    
    return matches, incorrect_arrow,step_correctness

def check_resonance(model, student, steps):
    # If no keys in 'student' end with '_a', return 'no_resonance'
    if not any(key.endswith('_ra') for key in student):
        return ['no_resonance']

    resonance = []
    for key in student:
        if key.endswith('_rb'):
            i = key.split('_')
            step_index = float(i[1])
            for x in steps['matching_sequences']:
                if (x[0] - x[2]) <= step_index <= (x[1] + 1 - x[2]):
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

def individual_steps(model,student,res_exe,reaction_arrows):
    
    ''' 
    The student answer of the reactin mechanism is assessed on a global level.
    First, individual reaction steps from both the model answer and student
    answer are generated as reaction transformations in SMILES. Next, the 
    individual steps between the model and student answer are compared. If 
    all steps between the student and model match, the mechanism is correct
    on a global level. Otherwise the inconsistencies are returned in a 
    dictionary.
    
    '''

    model_reactions = generate_reactions(model, reaction_type = 'model')
    student_reactions = generate_reactions(student, reaction_type = 'student')

    default = 'model_default'
    if len(student_reactions[0]) < len(model_reactions[0]):
        default = 'student_default'
    steps = compare_reactions(model_reactions[0],student_reactions[0])

    arrow_check = reaction_arrow_check(steps['individual_steps'],reaction_arrows)
    
    s_keys = steps['individual_steps'].keys()
    
    steps_bool = []
    for x in s_keys:
        steps_bool.append(steps['individual_steps'][x][0])

    exersice = True
    if all(all(not item for item in sublist[0]) for sublist in steps_bool): 
        model_m = []
        try:
            for x in list(model['reactants']['molecules'].keys()):
                model_m.append(model['reactants']['molecules'][x].replace('[H]',''))
        except:
            for x in list(model['reactants_ma']['molecules'].keys()):
                model_m.append(model['reactants_ma']['molecules'][x].replace('[H]',''))
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
        if res_exe == 'yes' or res_exe == 'Yes' or res_exe == 'YES':
            presence_resonance = check_resonance(model,student,steps)
            resonance = True
        else: 
            resonance = False

    if steps['index_resonance'] is False:
        return steps,model_reactions[0],student_reactions[0],exersice,list(model.keys()),{'resonance':presence_resonance,'resonance_present':resonance, 'default':default, 'reaction_arrows':arrow_check}
    else:
        return steps,model_reactions[0][steps['index_resonance']],student_reactions[0],exersice,model_reactions[1][steps['index_resonance']],{'resonance':presence_resonance,'resonance_present':resonance, "default":default, 'reaction_arrows':arrow_check}

def reaction_transformations(model,student,steps):

    if steps[-1]['default'] == 'student_main':
        short = student
        long = model
        short_keys = student_keys = sorted([item for item in list(student.keys()) if not item.endswith(("_rb"))], key=custom_reaction_sort)
        long_keys = model_keys = sorted(steps[4], key=custom_reaction_sort)
    else:
        short = model
        long = student
        short_keys = model_keys = sorted(steps[4], key=custom_reaction_sort)
        long_keys = student_keys = sorted([item for item in list(student.keys()) if not item.endswith(("_rb"))], key=custom_reaction_sort)
        
    seq = steps[0]['matching_sequences']
    transformations = []
    incorrect_arrow = {}
    step_correct = {}
    first_step = True

    len_steps = len(long_keys) - 1

    if seq[0][0] != 0:
        m = short[short_keys[0]].get('arrows')
        s = long['reactants'].get('arrows')
        check = compare_arrows(m,s,steps[1][0],steps[2][0])
        transformations.append(check)
        first_step = False
    
    for series in seq:
        for step in np.arange(series[0],series[1]+2, 1):
            if step == 0 and first_step == True:
                m = short[short_keys[0]]['arrows']
                if series[2] == 0:
                    s = long['reactants']['arrows']
                else:
                    comp = long_keys[step - series[2]]
                    s = long[comp]['arrows']

                check = compare_arrows(m,s,steps[1][step-series[2]],steps[2][step])
                transformations.append(check[0])
                incorrect_arrow[step] = check[1]
                step_correct[step] = check[2]
                if len(list(model.keys())) < 3:
                    break
            
            elif step > 0 and step < len_steps: ### should also be included if step == len_steps? Then reconsider how len_steps is defined depending on which mechanism is longest
                comp_m = short_keys[step-(series[2])] 
                if series[2] == 0:
                    comp_s = long_keys[step]
                elif series[2] != 0:
                    comp_s = long_keys[step]

                m = short[comp_m]['arrows']
                s = long[comp_s]['arrows']

                try:
                    check = compare_arrows(m,s,steps[1][step-series[2]],steps[2][step])
                except: 
                    if m:
                        check = compare_arrows(m,s,">>",">>")

                transformations.append(check[0])
                incorrect_arrow[step] = check[1]
                step_correct[step] = check[2]
            elif step >= len_steps:
                break
    return transformations, incorrect_arrow, step_correct
            
