from exersices_ORC import exersices
from convert_json_to_smiles import parse_json_to_smiles
from mechanism_analysis import individual_steps, reaction_transformations
from primary_feedback import analysis_feedback
from openai_input import generate_feedback

def analyze(openai_key:False,exersice,student_json_code):
    ''' This function operates the analysis from the ORC_reaction_analysis input '''

    ''' retrieve model answer '''    
    json_code = exersices(exersice)
    model = parse_json_to_smiles(json_code)
    
    ''' compile student answer '''
    student = parse_json_to_smiles(student_json_code,True)
    mol_struc = student[1]
    student = student[0]
    ''' global reaction mechanism comparison '''
    reaction_steps = individual_steps(model,student)

    ''' detailed reaction mechainsm comparison '''
    if len(reaction_steps[0]['matching_sequences']) > 0 and reaction_steps[3] == True:
        transformations = reaction_transformations(model,student,reaction_steps)
    else:
        transformations = []
    ''' primary analysis '''
    issues = analysis_feedback(reaction_steps,transformations,mol_struc)
    
    ''' formalize student feedback from openai analysis '''
    feedback = generate_feedback(openai_key,model,student,reaction_steps,transformations,mol_struc,issues)

    return feedback