import numpy as np
from rdkit import Chem
from .mechanism_comparison import compare_reactions
from .helper_functions import unique_strings, common_substring, flatten, filter_variations_main_res, filter_variations_res, custom_reaction_sort, _split_reaction
from .molecular_structures import are_constitutional_isomers, get_molecule_name

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
                    
                    reactants_str = ".".join(m for m in reactant_set)
                    products_str = ".".join(m for m in product_set)
        
                    reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
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
                    
                    reactants_str = ".".join(m for m in reactant_set)
                    products_str = ".".join(m for m in product_set)
        
                    reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
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
                
                reactants_str = ".".join(m for m in reactant_set)
                products_str = ".".join(m for m in product_set)

                reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
                reaction_list.append(reaction)
            
        elif resonance == False:  
            for i, step in enumerate(steps[:-1]):
                reactants = mols[steps[i]]
                products = mols[steps[i+1]]

                # Identify unique molecules in reactants and products
                reactant_set = set(reactants)
                product_set = set(products)

                reactants_str = ".".join(m for m in reactant_set)
                products_str = ".".join(m for m in product_set)
                reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"

                reaction_list.append(reaction)
                
    elif n_steps == 2:
        reactants = mols[steps[0]]
        products = mols[steps[1]]

        # Identify unique molecules in reactants and products
        reactant_set = set(reactants)
        product_set = set(products)
        
        reactants_str = ".".join(m for m in reactant_set)
        products_str = ".".join(m for m in product_set)

        reaction = f"{reactants_str}>>{products_str}" if products_str else f"{reactants_str}>>"
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
            if '[H-]' != a[2]:
                reaction_type.append('proton transfer')
                continue
            elif '[H-]' == a[2]:
                reaction_type.append('nucleophilic attack')
                continue
            
            ''' identify nucleophilic attack '''
        else:
            if a[4] == a[5]:
                if max(s.GetNumAtoms() for s in r_chem) > max(s.GetNumAtoms() for s in p_chem):
                    reaction_type.append('loss of leaving group')
                    continue
                else:
                    index_largest_reactant = max(enumerate(r_chem), key=lambda x: x[1].GetNumAtoms())[0]
                    index_largest_product = max(enumerate(p_chem), key=lambda x: x[1].GetNumAtoms())[0]
                    rearrange = are_constitutional_isomers(Chem.MolToSmiles(r_chem[index_largest_reactant]), Chem.MolToSmiles(p_chem[index_largest_product]))[0]
                    if rearrange == True:
                        reaction_type.append('rearrangement')
                    else:
                        reaction_type.append('internal electron movement')
                    continue

            elif a[4] != a[5]:
                reaction_type.append('nucleophilic attack')
                continue 

    return reaction_type

def reaction_arrow_check(steps,arrows):
    step_numbers = steps.keys()
    wrong_arrows = {}

    for s in step_numbers:
        if steps[s][0][0][0] == True:
            m_step = steps[s][0][1]
            if m_step == 'not_present':
                continue
            if arrows['s_arrows'][s] != arrows['m_arrows'][m_step] and arrows['m_arrows'][m_step] == 'synthetic':
                wrong_arrows[int(s)+1] = False
    return wrong_arrows

def compare_arrows(model_arrows, student_arrows,model_reaction,student_reaction,allowed_acid_base,all_acid_base,step_key):
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
            return ([reaction_rules(model_arr,model_reaction),[False,False]],False,False,False)   
        else:
            return ([[False, False],[False,False]],False,False,False)

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
    step_appearance = False
    matches = []
    max_it = len(short_arr) - 1
    common_acid = []
    incorrect_ab = {}

    for ind,x in enumerate(long_arr):
        app = False
        if student == 'short_mechanism' and len(short_arr) == len(matches):
            break 
        candidate_m = False
        for i,y in enumerate(short_arr):
            base_source_student = False
            base_source_model   = False
            base_source_flipped = False
            wrong_base_source_student = False
            acid_target_student = False
            acid_target_model   = False
            acid_target_flipped = False
            wrong_acid_target_student = False
            
            def m_obj(smiles):
                try:
                    m = Chem.MolFromSmiles(smiles)
                    return Chem.MolToSmiles(m)
                except:
                    return smiles
            
            if student == 'short_mechanism':
                for acid in allowed_acid_base[0]:
                    if m_obj(acid) == m_obj(y[5]):
                        acid_target_student = True
                    if m_obj(acid) == m_obj(x[5]):
                        acid_target_model   = True
                    if m_obj(acid) == m_obj(y[4]):
                        acid_target_flipped = True
                for base in allowed_acid_base[1]:
                    if m_obj(base) == m_obj(y[4]):
                        base_source_student = True
                    if m_obj(base) == m_obj(x[4]):
                        base_source_model = True
                    if m_obj(base) == m_obj(y[5]):
                        base_source_flipped = True
                for acid in all_acid_base[0]:
                    if m_obj(acid) == m_obj(y[5]):
                        wrong_acid_target_student = True
                for base in all_acid_base[1]:
                    if m_obj(base) == m_obj(y[4]):
                        wrong_base_source_student = True
            else:
                for acid in allowed_acid_base[0]:
                    if m_obj(acid) == m_obj(x[5]):
                        acid_target_student = True
                    if m_obj(acid) == m_obj(y[5]):
                        acid_target_model   = True
                    if m_obj(acid) == m_obj(x[4]):
                        acid_target_flipped = True
                for base in allowed_acid_base[1]:
                    if m_obj(base) == m_obj(x[4]):
                        base_source_student = True
                    if m_obj(base) == m_obj(y[4]):
                        base_source_model   = True
                    if m_obj(base) == m_obj(x[5]):
                        base_source_flipped = True
                for acid in allowed_acid_base[0]:
                    if m_obj(acid) == m_obj(x[5]):
                        wrong_acid_target_student = True
                for base in all_acid_base[1]:
                    if m_obj(base) == m_obj(x[4]):
                        wrong_base_source_student = True


            # correct arrow
            if [x[0],x[1],x[2],x[3]] == [y[0],y[1],y[2],y[3]]:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append((True,True,reaction_rules([x],model_reaction)[0]))
                else:
                    matches.append((True,True,reaction_rules([y],model_reaction)[0]))
                app = True
                break
            elif [x[1],x[3]] == [y[1],y[3]] and base_source_student == True and base_source_model == True:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append((True,True,reaction_rules([x],model_reaction)[0]))
                else:
                    matches.append((True,True,reaction_rules([y],model_reaction)[0]))
                app = True
                break
            elif [x[1],x[3]] == [y[1],y[3]] and wrong_base_source_student == True and base_source_model == True:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append((True,True,reaction_rules([x],model_reaction)[0]))

                    incorrect_ab[step_key] = ['incorrect base',get_molecule_name(y[2])]
                else:
                    matches.append((True,True,reaction_rules([y],model_reaction)[0]))
                    incorrect_ab[step_key] = ['incorrect base',get_molecule_name(x[2])]
                app = True
                break
            elif [x[0],x[2]] == [y[0],y[2]] and acid_target_student == True and acid_target_model == True:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append((True,True,reaction_rules([x],model_reaction)[0]))
                    if y[3] == "[H+]":
                        common_acid.append(True)
                else:
                    matches.append((True,True,reaction_rules([y],model_reaction)[0]))
                    if x[3] == "[H+]":
                        common_acid.append(True)
                app = True
                break
            elif [x[0],x[2]] == [y[0],y[2]] and wrong_acid_target_student == True and acid_target_model == True:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append((True,True,reaction_rules([x],model_reaction)[0]))
                    incorrect_ab[step_key] = ['incorrect acid',get_molecule_name(y[2])]
                else:
                    matches.append((True,True,reaction_rules([y],model_reaction)[0]))
                    incorrect_ab[step_key] = ['incorrect acid',get_molecule_name(x[2])]
            # correct start arrow
            elif (x[0],x[2]) == (y[0],y[2]):
                
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    candidate_m = (True,False,reaction_rules([x],model_reaction)[0])
                else:    
                    candidate_m = (True,False,reaction_rules([y],model_reaction)[0])
                app = True
                
            elif x[0] == y[0] and base_source_student == True and base_source_model == True:
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    candidate_m = (True,False,reaction_rules([x],model_reaction)[0])
                else:    
                    candidate_m = (True,False,reaction_rules([y],model_reaction)[0])
                app = True
                
            
            # correct end arrow
            elif (x[1],x[3]) == (y[1],y[3]):
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    candidate_m = (False,True,reaction_rules([x],model_reaction)[0])
                else:
                    candidate_m = (False,True,reaction_rules([y],model_reaction)[0])
                app = True
                
            elif x[1] == y[1] and acid_target_student == True and acid_target_model == True:
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    candidate_m = (False,True,reaction_rules([x],model_reaction)[0])
                else:    
                    candidate_m = (False,True,reaction_rules([y],model_reaction)[0])
                app = True
                
            
            # arrow flipped
            elif [x[0],x[1],x[2],x[3]] == [y[1],y[0],y[3],y[2]]:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append(('arrow_flipped',False,reaction_rules([x],model_reaction)[0]))
                else:
                    matches.append(('arrow_flipped',False,reaction_rules([y],model_reaction)[0]))
                app = True
                break

            elif [x[0],x[1],x[3]] == [y[1],y[0],y[2]] or [x[0],x[1],x[2]] == [y[1],y[0],y[3]] and acid_target_flipped == True and acid_target_model == True:
                step_appearance = True
                if x[6] != y[6]:
                    incorrect_arrow = True
                if student == 'short_mechanism':
                    matches.append(('arrow_flipped',False,reaction_rules([x],model_reaction)[0]))
                    if y[2] == "[H+]":
                        common_acid.append(True)
                else:
                    matches.append(('arrow_flipped',False,reaction_rules([y],model_reaction)[0]))
                    if x[2] == "[H+]":
                        common_acid.append(True)
                app = True
                break
            
            # arrow not present
            if app == False and max_it == i:
                if student == 'short_mechanism'  and len(long_mech) != len(short_mech):
                    continue
                else:
                    matches.append((False,False,'no_match'))
                    
            if app == True and max_it == i and candidate_m != False:
                matches.append(candidate_m)
                
    no_mistake = []
    for x in matches:
        no_mistake.append(x[0])
        no_mistake.append(x[1])
    step_correct = False
    
    if all(x == True for x in no_mistake) and len(long_arr) == len(short_arr):
        step_correct = True
        
    acid_values = {(True,True,'proton transfer'),(False,False,'proton transfer')}
    acid_values_flipped = {('arrow_flipped',False,'proton transfer'),(False,False,'proton transfer')}
    
    
    if True in common_acid and set(matches).issubset(acid_values):
        step_correct = True
        matches_and_model = [(True, True, 'proton transfer'),['proton transfer']]
    
    elif common_acid == True and set(matches).issubset(acid_values_flipped):
        matches_and_model = [("arrow_flipped", 'proton transfer'),['proton transfer']]
        
    elif student == 'short_mechanism':
        model_operation = reaction_rules(long_arr,model_reaction)
        matches_and_model = [matches,model_operation]
        
    elif student == 'long_mechanism':
        model_operation = reaction_rules(short_arr,model_reaction)
        matches_and_model = [matches,model_operation]

    return matches_and_model, incorrect_arrow,step_appearance,step_correct

def check_resonance(model, student, steps):
    # If no keys in 'student' end with '_a', return 'no_resonance'
    if not any(key.endswith('_ra') for key in student):
        return ['no_resonance']

    resonance = []
    for key in student:
        if key.endswith('_rb'):
            i = key.split('_')
            step_index = float(i[1])
            model_resonance = []
            for x in steps['matching_sequences']:
                if (x[0] - x[2]) <= step_index <= (x[1] + 1 - x[2]):
                    m_int = f"intermediates_{int(step_index) - x[2]}"
                    # Find all relevant model keys
                    m_res = [k for k in model if m_int in k]
                    # Extract molecules from model steps
                    step_molecules = [list(model[step]['molecules'].values()) for step in m_res]
                    # Find unique resonance molecules
                    resonance_molecules = unique_strings(step_molecules)
                    model_resonance.extend(resonance_molecules)

            # Get student's resonance molecules
            resonance_student = list(student[key]['molecules'].values())

            # Check if student's resonance molecule matches any valid ones
            if model_resonance and resonance_student[0] in model_resonance:
                resonance.append([True,key])
            else:
                resonance.append([False,key])
                
    return resonance

def individual_steps(model,student,res_exe,reaction_arrows,allowed_acid_base):
    
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
    steps = compare_reactions(model_reactions[0],student_reactions[0],allowed_acid_base)

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
            if x in student_m and x not in allowed_acid_base[0] and x not in allowed_acid_base[1]:
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
            presence_resonance = check_resonance(model,student,steps)

    if steps['index_resonance'] is False:
        return {'steps_check':steps,'reactions_model':steps['reactions_model'],'reactions_student':steps['reactions_student'],'correct_exersice':exersice,'model_keys':list(model.keys()),'secondary_check':{'resonance':presence_resonance,'resonance_present':resonance, 'default':default, 'reaction_arrows':arrow_check},'rogue_species':steps['rogue_species']}
    else:
        return {'steps_check':steps,'reactions_model':steps['reactions_model'],'reactions_student':steps['reactions_student'],'correct_exersice':exersice,'model_keys':model_reactions[1][steps['index_resonance']],'secondary_check':{'resonance':presence_resonance,'resonance_present':resonance, "default":default, 'reaction_arrows':arrow_check},'rogue_species':steps['rogue_species']}

def reaction_transformations(model,student,steps,allowed_acid_base,all_acid_base):

    if steps['secondary_check']['default'] == 'student_main':
        short = student
        long = model
        short_keys = student_keys = sorted([item for item in list(student.keys()) if not item.endswith(("_rb"))], key=custom_reaction_sort)
        long_keys = model_keys = sorted(steps['model_keys'], key=custom_reaction_sort)
    else:
        short = model
        long = student
        short_keys = model_keys = sorted(steps['model_keys'], key=custom_reaction_sort)
        long_keys = student_keys = sorted([item for item in list(student.keys()) if not item.endswith(("_rb"))], key=custom_reaction_sort)
   
    seq = steps['steps_check']['matching_sequences']
    transformations = {}
    incorrect_arrow = {}
    step_appearance = {}
    step_correct = {}
    
    first_step = True

    len_steps = len(long_keys) - 1

    if seq[0][0] != 0:
        m = short[short_keys[0]].get('arrows')
        s = long['reactants'].get('arrows')
        check = compare_arrows(m,s,steps['reactions_model'][0][0],steps['reactions_student'][0],allowed_acid_base,all_acid_base,0)
        transformations[0] = check
        first_step = False
    
    def generate_steps(series,l_model):
        final_step = series[1] - series[2]
        if final_step < l_model:
            list_steps = np.arange(series[0],series[1]+2, 1)
        else:
            list_steps = np.arange(series[0],series[1]+1, 1)
        return list_steps

    l_model = len(model)

    for iteration,series in enumerate(seq):
        for step in generate_steps(series,l_model):
            if step == 0 and first_step == True:
                m = short[short_keys[0]]['arrows']
                if series[2] == 0:
                    s = long['reactants']['arrows']
                else:
                    comp = long_keys[step - series[2]]
                    s = long[comp]['arrows']

                check = compare_arrows(m,s,steps['reactions_model'][0][step-series[2]],steps['reactions_student'][step],allowed_acid_base,all_acid_base,step)
                transformations[step]   = check[0]
                incorrect_arrow[step]   = check[1]
                step_appearance[step]   = check[2]
                step_correct[step]      = check[3]
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
                    check = compare_arrows(m,s,steps['reactions_model'][0][step-series[2]],steps['reactions_student'][step],allowed_acid_base,all_acid_base,step)
                except: 
                    if m:
                        check = compare_arrows(m,s,">>",">>",allowed_acid_base,all_acid_base,step)
                transformations[step]   = check[0]
                incorrect_arrow[step]   = check[1]
                step_appearance[step]   = check[2]
                step_correct[step]      = check[3]
            elif step >= len_steps:
                break

    return transformations, incorrect_arrow, step_appearance, step_correct
            
