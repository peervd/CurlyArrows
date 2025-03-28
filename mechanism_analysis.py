import numpy as np
from rdkit import Chem
from mechanism_comparison import compare_reactions, canonicalize_reaction

def generate_reactions(reaction_dict):
    reaction_list = []
    steps = list(reaction_dict.keys())
    n_steps = len(steps)

    # Extract molecules for each step
    mols = {s: list(reaction_dict[s]['molecules'].values()) for s in steps}

    if n_steps > 2:
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

    return reaction_list


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
        if 'H' in a[0][0] or 'H' in a[1][0]:
            reaction_type.append('proton transfer')
            continue

        ''' identify nucleophilic attack '''
        if 'H' not in a[0][0] and 'H' not in a[1][0] and a[2] != a[3]:
            reaction_type.append('nucleophilic attack')
        
        ''' identify loss of leaving group / rearrangemnt'''
        if 'H' not in a[0][0] and 'H' not in a[1][0] and a[2] == a[3]:

            if max(s.GetNumAtoms() for s in r_chem) > max(s.GetNumAtoms() for s in p_chem):
                reaction_type.append('loss of leaving group')
            elif max(s.GetNumAtoms() for s in r_chem) == max(s.GetNumAtoms() for s in p_chem):
                reaction_type.append('rearrangement')
            else:
                reaction_type.append('internal electron movement')
        
        #''' identify redox '''

    return reaction_type


def compare_arrows(model_arrows, student_arrows,model_reaction,student_reaction):
    if not model_arrows or not student_arrows:  # Handle empty input
        return [False, False]

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
        (long_mech[x].get('substructure_end'))]
        long_arr.append(arrows)
    short_arr = []
    for x in range(0,len(short_mech)):
        arrows = [(short_mech[x].get('start_arrow')),
        (short_mech[x].get('end_arrow')),
        (short_mech[x].get('substructure_start')),
        (short_mech[x].get('substructure_end'))]
        short_arr.append(arrows)
    match = []
    max_it = len(short_arr) - 1
    for x in long_arr:
        app = False
        for i,y in enumerate(short_arr):
            if x == y:
                match.append([ True,True ])
                app = True
                continue
            elif (x[0],x[2]) == (y[0],y[2]):
                match.append([ True,False ])
                app = True
                continue
            elif [x[1],x[3]] == [y[1],y[3]]:
                match.append([ False,True ])
                app = True
                continue
            elif x == [y[1],y[0],y[3],y[2]]:
                match.append([ 'arrow_flipped',False ])
                app = True
                continue
            elif app == False and max_it == i:
                match.append([ False,False ])


    if any(x == False for y in match for x in y):
        if student == 'short_mechanism':
            model_operation = reaction_rules(long_arr,model_reaction)
            student_operation = reaction_rules(short_arr,student_reaction)
            match = [model_operation,student_operation]


        elif student == 'long_mechanism':
            model_operation = reaction_rules(short_arr,model_reaction)
            student_operation = reaction_rules(long_arr,student_reaction)
            match = [model_operation,student_operation]

    return match


def individual_steps(model,student):
    
    ''' 
    The student answer of the reactin mechanism is assessed on a global level.
    First, individual reaction steps from both the model answer and student
    answer are generated as reaction transformations in SMILES. Next, the 
    individual steps between the model and student answer are compared. If 
    all steps between the student and model match, the mechanism is correct
    on a global level. Otherwise the inconsistencies are returned in a 
    dictionary.
    
    '''
    model_reactions = generate_reactions(model)
    student_reactions = generate_reactions(student)
    steps = compare_reactions(model_reactions,student_reactions)
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
            if x in student_m:
                correct_exercise.append(True)
            else:
                correct_exercise.append(False)

        if all(x == False for x in correct_exercise):
            exersice = False

    return steps,model_reactions,student_reactions,exersice

def reaction_transformations(model,student,steps):
    seq = steps[0]['matching_sequences']
    transformations = []
    len_steps = len(student) - 1
    first_step = False

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
                    comp = 'intermediates_%s' %(series[2]+1)
                    s = student[comp]['arrows']
                check = compare_arrows(m,s,steps[1][step],steps[2][step])
                
                if len(list(model.keys())) < 3:
                    break
            
            elif step > 0 and step < len_steps: ### if error: series[2] changed to series [1]
                comp_m = 'intermediates_%s' %(step-(series[2]))
                if series[2] == 0:
                    comp_s = 'intermediates_%s' %(step)
                elif series[2] != 0:
                    comp_s = 'intermediates_%s' %(step)
                m = model[comp_m]['arrows']
                s = student[comp_s]['arrows']

                check = compare_arrows(m,s,steps[1][step-series[2]],steps[2][step])
            elif step >= len_steps:
                break
            transformations.append(check)
  
    return transformations
            
