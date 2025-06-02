from exersices_ORC import exersices
from convert_json_to_smiles import parse_json_to_smiles
from mechanism_analysis import individual_steps, reaction_transformations
from primary_feedback import analysis_feedback
from openai_input import generate_feedback, generate_error_mechanism

def analyze(openai_key:False,exersice,student_json_code):
    ''' This function operates the analysis from the ORC_reaction_analysis input '''

    ''' retrieve model answer '''    
    exersice_info = exersices(exersice)
    code_key = list(exersice_info.keys())
    model = parse_json_to_smiles(exersice_info[code_key[0]],mech_category='model')
    model_reaction_arrow = model[1]
    model = model[0]
    
    ''' compile student answer '''
    student = parse_json_to_smiles(student_json_code,mech_category='student')

    if student[0] == 'invalid_mech':
        return generate_error_mechanism(student[1])
    
    reaction_arrows = {'m_arrows': model_reaction_arrow,'s_arrows': student[2]}
    
    mol_struc = student[1]
    student = student[0]
    
    ''' global reaction mechanism comparison '''
    reaction_steps = individual_steps(model,student,exersice_info.get('resonance'),reaction_arrows)

    ''' detailed reaction mechainsm comparison '''
    if len(reaction_steps[0]['matching_sequences']) > 0 and reaction_steps[3] == True:
        transformations = reaction_transformations(model,student,reaction_steps)
    else:
        transformations = []
    ''' primary analysis '''
    issues = analysis_feedback(reaction_steps,transformations,mol_struc)
    ''' formalize student feedback from openai analysis '''
    return generate_feedback(openai_key,model,student,reaction_steps,transformations,mol_struc,issues)
