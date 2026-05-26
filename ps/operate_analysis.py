from .exersices_ORC import exersices, get_acid_base, get_concept_tags, acid_base
from .convert_json_to_smiles import parse_json_to_smiles
from .mechanism_analysis import individual_steps, reaction_transformations
from .primary_feedback import analysis_feedback
from .openai_input import  generate_feedback, generate_error_mechanism, create_concept_tags, generate_prompt_rephrase, generate_prompt_categorize, generate_prompt_catC, generate_prompt_catD, assign_category
from .openai_communication import communicate_prompt
import os

def analyze(openai_key:False,exersice,student_json_code,student_reasoning:False):
    ''' This function operates the analysis from the ORC_reaction_analysis input '''

    ''' retrieve model answer '''    
    exersice_info = exersices(exersice)
    code_key = list(exersice_info.keys())
    model = parse_json_to_smiles(exersice_info[code_key[0]],mech_category='model')
    model_reaction_arrow = model[1]
    model = model[0]
    
    ''' acids and bases '''
    allowed_acid_base = get_acid_base(model,'list_acid_base.csv')
    all_acid_base = acid_base('list_acid_base.csv')

    ''' compile student answer '''
    student = parse_json_to_smiles(student_json_code,mech_category='student')

    if student[0] == 'invalid_mech':
        return generate_error_mechanism(student[1])
    
    reaction_arrows = {'m_arrows': model_reaction_arrow,'s_arrows': student[2]}
    
    mol_struc = student[1]
    student = student[0]

    ''' global reaction mechanism comparison '''
    reaction_steps = individual_steps(model,student,exersice_info.get('resonance'),reaction_arrows,allowed_acid_base)

    ''' detailed reaction mechainsm comparison '''
    if len(reaction_steps['steps_check']['matching_sequences']) > 0 and reaction_steps['correct_exersice'] == True:
        transformations = reaction_transformations(model,student,reaction_steps,allowed_acid_base,all_acid_base)
    else:
        transformations = []
    ''' primary analysis '''
    issues = analysis_feedback(reaction_steps,transformations,mol_struc)
    ''' formalize student feedback from openai analysis '''
    if student_reasoning and len(str(student_reasoning)) < 5:
        student_reasoning = False

    if openai_key == False :
        return generate_feedback(model,student,reaction_steps,transformations,mol_struc,issues)['feedback']
    elif student_reasoning == False:
        concept_tags_text = create_concept_tags(get_concept_tags('reaction_concept_tags.csv',exersice))
        feedback_CA = generate_feedback(model,student,reaction_steps,transformations,mol_struc,issues)
        if feedback_CA['category'] == 'A':
            return feedback_CA['feedback'] + "\n\nWrite down your reasoning for better feedback."
        else:
            prompt = generate_prompt_rephrase('ps/prompt_template_rephrase_CA.txt',feedback_CA['feedback'])
            return communicate_prompt(prompt,openai_key) + "\n\nWrite down your reasoning for better feedback."
    else:
        feedback_CA = generate_feedback(model, student, reaction_steps, transformations, mol_struc, issues)
        if feedback_CA['feed_AI'] == True:
            concept_tags_text = create_concept_tags(get_concept_tags('reaction_concept_tags.csv',exersice))
            prompt = generate_prompt_categorize('ps/prompt_template_categorize.txt',concept_tags_text,feedback_CA['feedback'],student_reasoning)
            primary_category = communicate_prompt(prompt,openai_key)
            #return primary_category #dont forget to remove feedback request!
            try:
                cat, letter = primary_category.split()
                feedback_category = assign_category(letter,feedback_CA['category'])
                if feedback_category == 'A':
                    return "Both your reaction mechanism and explanation are correct! Keep up the good work.\n\n" \
                    "This response was generated with the help of AI."
                elif feedback_category == 'B':
                    feedback = 'Your explanation for the reaction mechanism is correct. However, your reaction mechamism needs adjusting.'
                    feedback += feedback_CA['feedback']
                    return feedback
                elif feedback_category == 'C':
                    feedback_prompt = generate_prompt_catC('ps/prompt_template_catC.txt',concept_tags_text,student_reasoning)
                    return communicate_prompt(feedback_prompt,openai_key)
                elif feedback_category == 'D' or feedback_category == None:
                    feedback_prompt = generate_prompt_catD('ps/prompt_template_catD.txt',concept_tags_text,feedback_CA['feedback'],student_reasoning)
                    return communicate_prompt(feedback_prompt,openai_key)
            except:
                feedback_prompt = generate_prompt_catD('ps/prompt_template_general.txt',concept_tags_text,feedback_CA['feedback'],student_reasoning)
                return communicate_prompt(feedback_prompt,openai_key)
        else:
            return feedback_CA['feedback']
        