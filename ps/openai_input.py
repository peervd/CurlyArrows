from .helper_functions import list_to_string
from .molecular_structures import get_molecule_name
import re


def create_concept_tags(concept_tags):
    tag_prompt = "Here, we provide predefined background information on the reaction. "
    tag_prompt += "The reactants are %s and it produces %s. " %(concept_tags['reactants'],concept_tags['product'])
    tag_prompt += "It is a %s type of reaction with a total number of %s steps. " %(concept_tags['reaction_type'],len(concept_tags['steps'].keys()))
    tag_prompt += "Each step of the mechanism has a certain category of transformation and a number of concept tags. These are specified below: "
    for step in concept_tags['steps'].keys():
        tag_prompt += "Step %s of the mechanism is a %s. " %(step,concept_tags['steps'][step]['category'])
        tag_prompt += "The 'key' concept tags for step %s are: %s. 'Secondary' concept tags for step %s are: %s. " %(step,concept_tags['steps'][step]['tags_1'],step,concept_tags['steps'][step]['tags_2'])
    return tag_prompt
    
def generate_prompt_categorize(template_path,concept_tags,feedback_CA,student_reasoning):
    """
    Loads a markdown prompt template and fills in the three input sections.

    Parameters:
    - template_path (str): Path to the prompt_template.txt file.
    - concept_tags (str): Reaction Mechanism Overview (Reference).
    - feedback_CA (str): Student's Drawn Mechanism – Mechanical Analysis.
    - student_reasoning (str): Student’s Reasoning (Written Explanation).

    Returns:
    - str: The final filled-in prompt ready to be passed to the LLM.
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        template = f.read()
    # Replace the placeholders with the actual inputs
    template = template.replace(
        "**Paste instructor-provided mechanism summary here.**", concept_tags.strip())
    template = template.replace(
        "**Paste student’s explanation here.**", student_reasoning.strip())

    return template

def generate_prompt_rephrase(template_path,feedback_CA):
    """
    Loads a markdown prompt template and fills in the three input sections.

    Parameters:
    - template_path (str): Path to the prompt_template.txt file.
    - concept_tags (str): Reaction Mechanism Overview (Reference).
    - feedback_CA (str): Student's Drawn Mechanism – Mechanical Analysis.

    Returns:
    - str: The final filled-in prompt ready to be passed to the LLM.
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        template = f.read()
    # Replace the placeholders with the actual inputs
    template = template.replace(
        "*Paste mechanical analysis of student's drawing here.*", feedback_CA.strip())

    return template

def generate_prompt_catC(template_path,concept_tags,student_reasoning):
    """
    Loads a markdown prompt template and fills in the three input sections.

    Parameters:
    - template_path (str): Path to the prompt_template.txt file.
    - concept_tags (str): Reaction Mechanism Overview (Reference).
    - feedback_CA (str): Student's Drawn Mechanism – Mechanical Analysis.
    - student_reasoning (str): Student’s Reasoning (Written Explanation).

    Returns:
    - str: The final filled-in prompt ready to be passed to the LLM.
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        template = f.read()
    # Replace the placeholders with the actual inputs
    template = template.replace(
        "**Paste instructor-provided mechanism summary here.**", concept_tags.strip())
    template = template.replace(
        "**Paste student’s explanation here.**", student_reasoning.strip())

    return template

def generate_prompt_catD(template_path,concept_tags,feedback_CA,student_reasoning):
    """
    Loads a markdown prompt template and fills in the three input sections.

    Parameters:
    - template_path (str): Path to the prompt_template.txt file.
    - concept_tags (str): Reaction Mechanism Overview (Reference).
    - feedback_CA (str): Student's Drawn Mechanism – Mechanical Analysis.
    - student_reasoning (str): Student’s Reasoning (Written Explanation).

    Returns:
    - str: The final filled-in prompt ready to be passed to the LLM.
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        template = f.read()
    # Replace the placeholders with the actual inputs
    template = template.replace(
        "**Paste instructor-provided mechanism summary here.**", concept_tags.strip())
    template = template.replace(
        "**Paste mechanical analysis of student's drawing here.**", feedback_CA.strip())
    template = template.replace(
        "**Paste student’s explanation here.**", student_reasoning.strip())

    return template

def assign_category(AI_cat,CA_cat):
    if AI_cat == 'A' and CA_cat == 'A':
        return 'A'
    elif AI_cat == 'A' and CA_cat == 'B':
        return 'B'
    elif AI_cat == 'B' and CA_cat == 'A':
        return 'C'
    elif AI_cat == 'B' and CA_cat == 'B':
        return 'D'
    else:
        return None


def generate_feedback(model,student,reaction_steps,transformations,mol_struc,issues):
    AI = False
    if reaction_steps['correct_exersice'] == False:
        feedback = "Make sure you select the correct exercise. You do not start your mechanism with the correct reactants. "
        return {'feedback':feedback,'feed_AI':AI,'category':None}

    if all(reaction_steps['steps_check']['individual_steps'][x] == [([False, False], 'not_present')] for x in reaction_steps['steps_check']['individual_steps'].keys()) and not transformations:
        feedback = "None of the steps and curly arrows match with the model answer. Check if you did select the correct exercise. "
        return {'feedback':feedback,'feed_AI':AI,'category':None}
    AI = True
    """Generates structured feedback based on analysis."""
    single_feedback = False
    feedback = "Your input was carefully analyzed. "
    issue_list = list(issues.values())
    if all(x == True for x in issue_list):
        AI = True
        feedback = "Your reaction mechanism is correct. "
        return {'feedback':feedback,'feed_AI':AI,'category':'A'}
    
    elif issues["global"] == False:
        incorrect_steps = []
        s_keys = reaction_steps['steps_check']['individual_steps'].keys()
        steps_bool = []
        for x in s_keys:
            steps_bool.append(reaction_steps['steps_check']['individual_steps'][x][0])

        for i,x in enumerate(steps_bool):
            if x[0][0] == False and i == 0:
                incorrect_steps.append(i + 1)
            if x[0][1] == False:
                if (i + 1) not in incorrect_steps:
                    incorrect_steps.append(i + 1)
        if all(x == True for x in reaction_steps['steps_check']['product'].values()):
            feedback += "Although you do come with the correct product, the reaction mechanism does not follow the correct steps: "
        elif reaction_steps['steps_check']['product']['produced'] == True and reaction_steps['steps_check']['product']['last_step'] == False: 
            feedback += "Although you do produce the correct product, you incorrectly let it react further. "
        elif reaction_steps['steps_check']['product']['produced'] == False and reaction_steps['steps_check']['all_steps_match'] == False:
            feedback += "The reaction mechanism does not follow the correct steps and the correct product is not formed: "
        elif reaction_steps['steps_check']['product']['produced'] == False and reaction_steps['steps_check']['all_steps_match'] == True:
            feedback += "The reaction mechanism follows the correct steps, but the correct product is not formed: "    
        if len(reaction_steps['reactions_model'][0]) != len(reaction_steps['reactions_student']):
            feedback += "The model answer has %s intermediate steps whereas you draw %s . " %(len(reaction_steps['reactions_model'][0])+1, len(reaction_steps['reactions_student'])+1)
        combine_steps = False
        all_split_steps = []
        if '>>' in reaction_steps['reactions_student']: 
            for index,r in enumerate(reaction_steps['reactions_student']):
                if r == '>>':
                    split_steps = [index+1,index+2]
                    all_split_steps.append(index+1)
                    all_split_steps.append(index+2)
                    if combine_steps == True:
                        feedback += "Also, in your mechanism steps %s should probably be combined. " % list_to_string(split_steps)
                    else:
                        feedback += "In your mechanism steps %s should probably be combined. " % list_to_string(split_steps)
                    combine_steps = True
        if all_split_steps:
            incorrect_steps = [x for x in incorrect_steps if (x + 1)  not in all_split_steps]
        if len(incorrect_steps) == 1:
            if incorrect_steps[0] in reaction_steps['steps_check']['missing_steps']:
                pass
            else:
                i_step = [incorrect_steps[0]+1]
                single_feedback = i_step
                if len(s_keys)+1 == i_step[0]:
                    feedback += "Reconsider how you draw the product. "
                else:    
                    feedback += "Step %s includes incorrectly drawn molecule(s). " % list_to_string(i_step)
        elif len(incorrect_steps) > 1:
            if incorrect_steps in reaction_steps['steps_check']['missing_steps']:
                new_incorrect_steps = []
                for x in incorrect_steps:
                    if (x-1) in reaction_steps['steps_check']['missing_steps']: #maybe not x-1
                        continue
                    new_incorrect_steps.append(x+1)
                
                if len(new_incorrect_steps) == 1:
                    feedback += "Step %s includes incorrectly drawn molecule(s). " % list_to_string(new_incorrect_steps)
                elif len(new_incorrect_steps) > 1:
                    feedback += "Steps %s include incorrectly drawn molecule(s). " % list_to_string(new_incorrect_steps)
            else:    
                new_incorrect_steps = []
                for x in incorrect_steps:
                    new_incorrect_steps.append(x+1)
                feedback += "Steps %s in your drawing do not correspond with the model answer. " % list_to_string(new_incorrect_steps)

        if reaction_steps['steps_check']['missing_steps']:
            r_missing = []
            for x in reaction_steps['steps_check']['missing_steps']:
                if transformations[2].get(x-1) == False and x not in r_missing:
                    r_missing.append(str(x))
                r_missing.append(str(x+1))
            if len(r_missing) == 1:
                if r_missing[0] == '1':
                    feedback += "The first step from your mechanism is missing. "
                else:
                    feedback += "Step %s from the model answer does not appear in your mechanism. Reconsider how you continue after step %s." % (list_to_string(r_missing),str(reaction_steps['steps_check']['missing_steps'][0]))
                
            elif len(r_missing) > 1:
                if r_missing[0] == '1':
                    feedback += "The first step from your mechanism is missing. Also steps %s from the model answer do not appear in your mechanism. " % (list_to_string(r_missing.pop(0)))
                else:
                    feedback += "Steps %s from the model answer do no appear in your mechanism. Reconsider how you continue after step %s." % (list_to_string(r_missing),str(reaction_steps['steps_check']['missing_steps'][0]))
        if any(x == False for x in reaction_steps['secondary_check']['reaction_arrows'].values()):
            reaction_arrows = list(reaction_steps['secondary_check']['reaction_arrows'].keys())
            if len(reaction_arrows) == 1:
                feedback += "For step %s you should consider whether the step is reversible. " % list_to_string(reaction_arrows)
            elif len(reaction_arrows) > 1:
                feedback += "For steps %s you should consider whether the step is reversible. " % list_to_string(reaction_arrows)    

    if issues["mechanistic"] == False and issues["structure"] == False:

        same_step =  []
        for x in transformations[3]:
            if x in reaction_steps['steps_check']["molecular_structures"].keys():
                same_step.append(x+1)
        
        if same_step == single_feedback:
            feedback += "The molecular structure is incorrect and it is not possible to assess the electron movement. First make sure you draw the molecules correctly. " 
            step = same_step[0]-1
            for x in reaction_steps['steps_check']["molecular_structures"][step]:
                if transformations[-1][step-1] == False:
                    pass
                else:
                    feedback += "The %s carbon %s. "%(x[1],x[0])
        
        else:
            if same_step and len(same_step) == 1:
                feedback += "The molecular structure is step %s is incorrect and it is not possible to assess the electron movement. First make sure you draw the molecules correctly. " % list_to_string(same_step)
                step = same_step[0]-1
                for x in reaction_steps['steps_check']["molecular_structures"][step]:
                    if transformations[-1][step-1] == False:
                        pass
                    else:
                        feedback += "The %s carbon %s. "%(x[1],x[0])
            elif same_step:
                feedback += "The molecular structure in steps %s is incorrect and it is not possible to assess the electron movement. First make sure you draw the molecules correctly. " % list_to_string(same_step)
                step = same_step[0]-1
                for x in reaction_steps['steps_check']["molecular_structures"][step]:
                    if transformations[-1][step-1] == False:
                        pass
                    else:
                        feedback += "The %s carbon %s. "%(x[1],x[0])
            
        for x in same_step:
            reaction_steps['steps_check']["molecular_structures"].pop(x-1)
            transformations[3][x-1] = True

            

    if issues["mechanistic"] == False:
        if all(x == True for x in transformations[3].values()) and all(x == False for x in transformations[1].values()):
            pass
        else:
            try:
                for x in split_steps:
                    if transformations[2][x-1] == False and transformations[3][x-1] == False:
                        transformations[2].pop(x-1)
                        transformations[3].pop(x-1)
            except:
                pass
            report_arrows = False
            for step in transformations[3].keys():
                transformations[3][step]
                if transformations[3][step] == True:
                    continue
                elif step not in reaction_steps['steps_check']['missing_steps']:
                    report_arrows = True

            if report_arrows == True:                
                feedback += "Please pay attention to the electron flow (curly arrows):/n"

                for arrow in transformations[3].keys():
                    if transformations[3][arrow] == False:
                        mistakes_arrow = transformations[0][arrow][0]
                        model_arrows = transformations[0][arrow][1]
                        report_mistake = {'correct':0,'start correct':0,'does not start correct':0,'flipped':0,'not match':0}
                        correct_arrows = {}
                        false_transformations = []
                        for mistake in mistakes_arrow:
                            false_transformations.append(mistake[2])
                            if mistake[0] == True and mistake[1] == True:
                                report_mistake['correct'] +=1
                                continue
                            elif mistake[0] == 'arrow_flipped' and mistake[1] == False:
                                report_mistake['flipped'] +=1
                                continue
                            elif mistake[0] == True and mistake[1] == False:
                                report_mistake['start correct'] +=1
                                continue
                            elif mistake[0] == False and mistake[1] == True:
                                report_mistake['does not start correct'] +=1
                                continue
                            elif mistake[0] == False and mistake[1] == False:
                                report_mistake['not match'] +=1       
                
                        if 'internal electron movement' in model_arrows:
                            n_arrows = len(model_arrows)
                            
                            for trans in model_arrows:
                                if trans != "internal electron movement":
                                    general_arrow = trans
                            correct_arrows[arrow+1] = [general_arrow] * n_arrows
    
                        else:
                            correct_arrows[arrow+1] = model_arrows
    
                        ### report arrows
                        if sum(list(report_mistake.values())) == 0:
                            feedback += "Include curly arrows in step %s. " % (arrow+1)
                        
                        elif len(correct_arrows[arrow+1]) == report_mistake['correct']:
                            if report_mistake['correct'] == 1:
                                if sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]) == 1:
                                    feedback += "In step %s you draw the correct arrow. However, also %s arrow is incorrect. " % (arrow+1,sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]))
                                elif sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]) > 1:
                                    feedback += "In step %s you draw the correct arrow. However, also %s arrows are incorrect. " % (arrow+1,sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]))
                            else:  
                                if sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]) == 1:
                                    feedback += "In step %s you draw the correct arrows. However, also %s arrow is incorrect. " % (arrow+1,sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]))
                                elif sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]) > 1:
                                    feedback += "In step %s you draw correct arrows. However, also %s arrows are incorrect. " % (arrow+1,sum([report_mistake['not match'],report_mistake['start correct'],report_mistake['does not start correct'],report_mistake['flipped']]))
                            
                        elif sum(list(report_mistake.values())) == len(correct_arrows[arrow+1]):
                            if report_mistake['correct'] + 1 == len(correct_arrows[arrow+1]):
                                if set(false_transformations) == set(model_arrows):
                                    feedback += "Step %s contains curly arrows for the correct type of transformation, but the transformation is not correct. One incorrect arrow " %(arrow+1)
                                    for x in report_mistake.keys():
                                        if x == 'correct':
                                            continue
                                        if report_mistake[x] > 0:
                                            if x == 'flipped':
                                                feedback += 'is pointed in the wrong direction. '
                                            else:
                                                feedback += 'does %s. ' %x
                                else:
                                    feedback += "Step %s should include %s arrows for a %s. The one incorrect arrow " %(arrow+1,len(correct_arrows[arrow+1]),correct_arrows[arrow+1][0].replace('[','').replace(']',''))
                                    for x in report_mistake.keys():
                                        if x == 'correct':
                                            continue
                                        if report_mistake[x] > 0:
                                            if x == 'flipped':
                                                feedback += 'is pointed in the wrong direction. '
                                            else:
                                                feedback += 'does %s. ' %x
                            else:
                                if set(false_transformations) == set(model_arrows):
                                    feedback += "Step %s contains arrows for the correct type of transformations, but the transformation is not correct. Here, " %(arrow+1)
                                    for x in report_mistake.keys():
                                        if x == 'correct':
                                            continue
                                        if report_mistake[x] > 0:
                                            feedback += ', '
                                            if x == 'flipped':
                                                if report_mistake[x] == 1:
                                                    feedback += 'one incorrect arrow is pointed in the wrong direction '
                                                else:
                                                    feedback += '%s incorrect arrows are pointed in the wrong direction ' %report_mistake[x]
                                            else:
                                                if report_mistake[x] == 1:
                                                    feedback += 'one incorrect arrow does %s ' %x
                                                else:
                                                    feedback += '%s incorrect arrows do %s ' %(report_mistake[x],x)
                                    feedback += '. '
                                    
                                else:    
                                    feedback += "Step %s should include %s arrows for a %s. However " %(arrow+1,len(correct_arrows[arrow+1]),correct_arrows[arrow+1][0].replace('[','').replace(']',''))
                                    for x in report_mistake.keys():
                                        if x == 'correct':
                                            continue
                                        if report_mistake[x] > 0:
                                            feedback += ', '
                                            if x == 'flipped':
                                                if report_mistake[x] == 1:
                                                    feedback += 'one incorrect arrow is pointed in the wrong direction '
                                                else:
                                                    feedback += '%s incorrect arrows are pointed in the wrong direction ' %report_mistake[x]
                                            else:
                                                if report_mistake[x] == 1:
                                                    feedback += 'one incorrect arrow does %s ' %x
                                                else:
                                                    feedback += '%s incorrect arrows do %s ' %(report_mistake[x],x)
                                    feedback += '. '
                                          
                        elif sum(list(report_mistake.values())) < len(correct_arrows[arrow+1]):
                            
                            if sum(list(report_mistake.values())) == 1:
                                feedback += "Step %s should include %s arrow for a %s. However" %(arrow+1,len(correct_arrows[arrow+1]),correct_arrows[arrow+1][0].replace('[','').replace(']',''))
                                for x in report_mistake.keys():
                                    if x == 'correct' and report_mistake[x] > 0:
                                        feedback += ', you draw only one arrow correctly'
                                    elif report_mistake[x] > 0:
                                        feedback += ', '
                                        if x == 'flipped':
                                            feedback += 'you include only one incorrect arrow which is pointed in the wrong direction'
                                        else:
                                            feedback += 'you include only one incorrect arrow that does %s' %x                                
                                feedback += '. '
                                
                            elif sum(list(report_mistake.values())) > 1:
                                feedback += "Step %s should include %s arrows for a %s. However, you include only %s arrows. Here" %(arrow+1,len(correct_arrows[arrow+1]),correct_arrows[arrow+1][0].replace('[','').replace(']',''),sum(list(report_mistake.values())))
                                for x in report_mistake.keys():
                                    if x == 'correct':
                                        if report_mistake[x] == 1:
                                            feedback += ', one arrow is drawn correctly'
                                        else:
                                            feedback += ', %s arrows are drawn correctly' %report_mistake[x]
                                    elif report_mistake[x] > 0:
                                        feedback += ', '
                                        if x == 'flipped':
                                            
                                            if report_mistake[x] == 1:
                                                feedback += 'one incorrect arrow is pointed in the wrong direction'
                                            else:
                                                feedback += '%s incorrect arrows are pointed in the wrong direction' %report_mistake[x]
                                        else:
                                            if report_mistake[x] == 1:
                                                feedback += 'one incorrect arrow does %s' %x
                                            else:
                                                feedback += '%s incorrect arrows do %s' %(report_mistake[x],x)
                                feedback += '. '
                            
                        elif sum(list(report_mistake.values())) > len(correct_arrows[arrow+1]):
                            feedback += "Step %s should include %s arrows for a %s. However, you include only %s arrows. Here" %(arrow+1,len(correct_arrows[arrow+1]),correct_arrows[arrow+1][0].replace('[','').replace(']',''),sum(list(report_mistake.values())))
                            for x in report_mistake.keys():
                                if x == 'correct':
                                    if report_mistake[x] == 1:
                                        feedback += ', one arrow is drawn correctly'
                                    else:
                                        feedback += ', %s arrows are drawn correctly' %report_mistake[x]
                                elif report_mistake[x] > 0:
                                    feedback += ', '
                                    if x == 'flipped':
                                        
                                        if report_mistake[x] == 1:
                                            feedback += 'one incorrect arrow is pointed in the wrong direction'
                                        else:
                                            feedback += '%s incorrect arrows are pointed in the wrong direction' %report_mistake[x]
                                    else:
                                        if report_mistake[x] == 1:
                                            feedback += 'one incorrect arrow does %s' %x
                                        else:
                                            feedback += '%s incorrect arrows do %s' %(report_mistake[x],x)
                            feedback += '. '
            all_false = all(not any(v) if isinstance(v, list) else not v for v in transformations[1].values())
            if all_false == False:    
                all_steps = transformations[1].keys()
                false_arrow = []
                for x in all_steps:
                    if transformations[1][x] == True:
                        false_arrow.append(x+1)
                if len(false_arrow) == 1:
                    feedback += "In step %s you should reconsider the number of electrons indicated in the electron push (half or full arrow head). " % list_to_string(false_arrow)
                elif len(false_arrow) >1:
                    feedback += "In steps %s you should reconsider the number of electrons indicated in the electron push (half or full arrow head). " % list_to_string(false_arrow)
                    
    if issues["structure"] == False:
        
        if reaction_steps['steps_check']["molecular_structures"] and all(x == False for x in reaction_steps['steps_check']["molecular_structures"].values()):
            pass
        else:
            student_steps = sorted(mol_struc.keys(), key=lambda x: (0 if 'reactant' in x else (2 if 'product' in x else 1)))
            structure_comparison = reaction_steps['steps_check']["molecular_structures"]
            for i,x in enumerate(student_steps):
                if any(v is False for v in mol_struc[x].values()) or (i) in structure_comparison:
                    if x == 'reactants':
                        if len(mol_struc[x].keys()) == 1:
                            if (i) in structure_comparison:
                                feedback += "Check the molecular structure of your reactants. "
                                for x in structure_comparison[i]:
                                    feedback += "The %s. " % x[0]
                            else:
                                feedback += "Check the molecular structure of your reactants for valency errors. "
                        else:
                            if (i) in structure_comparison:
                                feedback += "Check the molecular structure of your reactants. "
                                for x in structure_comparison[i]:
                                    feedback += "The %s. " % x[0]
                            else:
                                feedback += "Check the molecular structure of your reactants for valency errors. "
                    if x == 'products':
                        if len(mol_struc[x].keys()) == 1:
                            if (i) in structure_comparison:
                                feedback += "Check the molecular structure of your products. "
                                for x in structure_comparison[i]:
                                    feedback += "The %s. " % x[0]
                            else:
                                feedback += "Check the molecular structure of your products for valency errors. "
                        else:
                            feedback += "Check the molecular structure of your products for valency errors. "
                    if 'intermediate' in x:
                        step = int(re.search(r"\d+", x).group()) +1
                        if sum(1 for value in mol_struc[x].values() if value is False) == 1 and not i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates in step %s for valency errors, as one is not correct. " %str(step)
                        elif sum(1 for value in mol_struc[x].values() if value is False) == 1 and i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s for valency errors. Here, the %s. " %(structure_comparison[i][0][1],str(step),structure_comparison[i][0][0])
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and not i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates in step %s for valency errors, as multiple are not correct. " %str(step)
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates with %s carbons in step %s for valency errors. Here, the %s. " %(structure_comparison[i][0][1],str(step),structure_comparison[i][0][0])
                        elif all(x == True for x in mol_struc[x].values()) and structure_comparison[i]:
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s. Here, the %s. " %(structure_comparison[i][0][1],str(step),structure_comparison[i][0][0])
                        
    if issues["resonance"] == False:
        if reaction_steps['secondary_check']['resonance_present'] == False:
            feedback += "In this exercise resonance stabilization does not play a significant role. You should remove the resonance structure from your mechanism. "
        else:
            feedback += "Include a correct resonance structure. "
            correct_res = []
            for x in reaction_steps['model_keys']:
                m_res = x.split('_')
                if 'ra' == m_res[-1] or 'rb' == m_res[-1]:
                    step_res = int(m_res[-2]) +1
                    correct_res.append(step_res)
            incorrect_steps = []
            for key in reaction_steps['secondary_check']['resonance']:
                if key == 'no_resonance' and reaction_steps['secondary_check']['resonance_present'] == True:
                    feedback += "One of the intermediates is stabalized by resonance. Make sure to draw this structure below the respective molecule. "
                    break
                key = key[1].split('_') 
                incorrect_steps.append(int(key[1])+1)
            if len(incorrect_steps) == 1:
                feedback += "In step %s you do not include the correct resonance structure. It should be included in step %s. " % (list_to_string(incorrect_steps),list_to_string(correct_res))
            elif len(incorrect_steps) > 1:
                feedback += "In steps %s you do not include the correct resonance structure. It should be included in steps %s. " % (list_to_string(incorrect_steps),list_to_string(correct_res))

            if reaction_steps['secondary_check']['resonance_present'] == False:
                feedback += "Check in the exersice whether you should include a resonance structure. "
    
    if issues['rogue'] == False:
        for x in reaction_steps['rogue_species'].keys():
            feedback += "In step %s you draw %s which does not belong here. "%(x+1,get_molecule_name(reaction_steps['rogue_species'][x][0]))
    feedback += "Good luck with the adjustments to your mechanism!"
    return {'feedback':feedback,'feed_AI':AI,'category':'B'}

def generate_error_mechanism(arrows):
    """ In case the mechanism does not follow the correct overall structure, a feedback is generated here."""

    wrong_arrows = []
    for a in arrows.keys():
        if arrows[a] == False:
            wrong_arrows.append(a)
    
    if len(wrong_arrows) == 1:
        feedback = "The mechanism you submitted does not follow the correct structure. You did not draw the %s arrow in the correct direction. Make sure you draw the mechanism horizontally, with horizontal reaction arrows separating each step. Resonance structures should be drawn below the main chain, with a vertical two headed resonance arrow. <br> Good luck with the adjustments to your mechanism! Paste the JSON code of your new reaction mechanism in the field above for a new analysis. " % list_to_string(wrong_arrows)
        return feedback      
    if len(wrong_arrows) > 1:
        feedback = "The mechanism you submitted does not follow the correct structure. You did not draw the %s arrows in the correct direction. Make sure you draw the mechanism horizontally, with horizontal reaction arrows separating each step. Resonance structures should be drawn below the main chain, with a vertical resonance arrow. <br> Good luck with the adjustments to your mechanism! Paste the JSON code of your new reaction mechanism in the field above for a new analysis. " % list_to_string(wrong_arrows)
        return feedback
