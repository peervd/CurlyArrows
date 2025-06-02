from helper_functions import list_to_string
import numpy as np
import openai
import requests
import re

def generate_prompt_chat_GPT(openai_key, model_reaction, student, reaction_steps, transformations, issues):
    # Initialize OpenAI client with API key
    client = openai.OpenAI(api_key=openai_key)

    # Define your prompt
    prompt_text = "Can you say - Hello world -?"

    # Call the API
    response = client.chat.completions.create(
        model='gpt-3.5-turbo',
        messages=[{"role": "user", "content": prompt_text}])
    return response.choices[0].message.content


def generate_prompt(hf_api_key, model, student, reaction_steps, transformations, mol_struc, issues):
    model_name = "facebook/opt-1.3b"

    API_URL = f"https://api-inference.huggingface.co/models/{model_name}"
    headers = {"Authorization": f"Bearer {hf_api_key}"}

    prompt_text = "Can you say - Hello world -?"
    payload = {"inputs": prompt_text, "parameters": {"max_new_tokens": 100}}

    response = requests.post(API_URL, headers=headers, json=payload)

    # Handle response errors
    if response.status_code == 200:
        return response.json()
    elif response.status_code == 403:
        return "Error 403: Invalid API token permissions. Check Hugging Face token settings."
    elif response.status_code == 404:
        return "Error 404: Model not found. Check model name on Hugging Face."
    else:
        return f"Error {response.status_code}: {response.text}"


def generate_feedback(openai_key,model,student,reaction_steps,transformations,mol_struc,issues):
    if reaction_steps[3] == False:
        return "Make sure you select the correct exercise. You do not start your mechanism with the correct reactants"
    
    """Generates structured feedback based on analysis."""
    if openai_key == False:
        feedback = "Your input was carefully analyzed and here is your feedback without openAI analysis.<br><br>"
        issue_list = list(issues.values())
        if all(x == True for x in issue_list):
            return "Your reaction mechanism is correct."
        
        elif issues["global"] == False:
            incorrect_steps = []
            s_keys = reaction_steps[0]['individual_steps'].keys()
            steps_bool = []
            for x in s_keys:
                steps_bool.append(reaction_steps[0]['individual_steps'][x][0])

            for i,x in enumerate(steps_bool):
                if x[0][0] == False and i == 0:
                    incorrect_steps.append(i + 1)
                if x[0][1] == False:
                    if (i + 1) not in incorrect_steps:
                        incorrect_steps.append(i + 1)
            if reaction_steps[0]['product'] == True:
                feedback += "Although you do come with the correct product, the reaction mechanism does not follow the correct steps:<br>"
            elif reaction_steps[0]['product'] == False:
                feedback += "The reaction mechanism does not follow the correct steps and the correct product is not formed:<br>"
            if len(reaction_steps[1]) != len(reaction_steps[2]):
                feedback += "The model answer has %s steps whereas you draw %s steps.<br>" %(len(reaction_steps[1])+1, len(reaction_steps[2])+1)
            combine_steps = False
            if '>>' in reaction_steps[2]:   
                for index,r in enumerate(reaction_steps[2]):
                    if r == '>>':
                        split_steps = [index+1,index+2]
                        for x in split_steps:
                            if x in incorrect_steps:
                                incorrect_steps.remove(x)
                                if combine_steps == True:
                                    feedback += "Also, in your mechanism steps %s should probably be combined.<br>" % list_to_string(split_steps)
                                feedback += "In your mechanism steps %s should probably be combined.<br>" % list_to_string(split_steps)
                                combine_steps = True
            if len(incorrect_steps) == 1:

                if incorrect_steps[0] in reaction_steps[0]['missing_steps']:
                    pass
                else:
                    i_step = [incorrect_steps[0]+1]
                    if len(s_keys)+1 == i_step[0]:
                        feedback += "Reconsider how you draw the product.<br>"
                    else:    
                        feedback += "Step %s in your drawing does not correspond with the model answer.<br>" % list_to_string(i_step)
            elif len(incorrect_steps) > 1:
                if incorrect_steps in reaction_steps[0]['missing_steps']:
                    new_incorrect_steps = []
                    for x in incorrect_steps:
                        if (x-1) in reaction_steps[0]['missing_steps']: #maybe not x-1
                            continue
                        new_incorrect_steps.append(x+1)
                    
                    if len(new_incorrect_steps) == 1:
                        feedback += "Step %s in your drawing does not correspond with the model answer.<br>" % list_to_string(new_incorrect_steps)
                    elif len(new_incorrect_steps) > 1:
                        feedback += "Steps %s in your drawing do not correspond with the model answer.<br>" % list_to_string(new_incorrect_steps)
                else:    
                    new_incorrect_steps = []
                    for x in incorrect_steps:
                        new_incorrect_steps.append(x+1)
                    feedback += "Steps %s in your drawing do not correspond with the model answer.<br>" % list_to_string(new_incorrect_steps)
            
            if reaction_steps[0]['missing_steps']:
                r_missing = []
                for x in reaction_steps[0]['missing_steps']:
                    if transformations[2].get(x-1) == False and x not in r_missing:
                        r_missing.append(str(x))
                    r_missing.append(str(x+1))
                if len(r_missing) == 1:
                    feedback += "Step %s from the model answer does not appear in your mechanism. Reconsider how you continue after step %s." % (list_to_string(r_missing),str(reaction_steps[0]['missing_steps'][0]))
                elif len(r_missing) > 1:
                    feedback += "Steps %s from the model answer do no appear in your mechanism. Reconsider how you continue after step %s." % (list_to_string(r_missing),str(reaction_steps[0]['missing_steps'][0]))
            if any(x == False for x in reaction_steps[-1]['reaction_arrows'].values()):
                reaction_arrows = list(reaction_steps[-1]['reaction_arrows'].keys())
                if len(reaction_arrows) == 1:
                    feedback += "For step %s you should consider whether the step is reversible." % list_to_string(reaction_arrows)
                elif len(reaction_arrows) > 1:
                    feedback += "For steps %s you should consider whether the step is reversible." % list_to_string(reaction_arrows)    
    
        if issues["mechanistic"] == False:
            report_arrows = False
            for i,arrows in enumerate(transformations[0]):
                if all(x is True for a in arrows for x in a):
                    continue
                elif i not in reaction_steps[0]['missing_steps']:
                    report_arrows = True
            if report_arrows == True:
                feedback += "<br>Please pay attention to the electron flow (curly arrows):<br>"
                if any(x == False or isinstance(x, str) for z in transformations[0] for y in z for x in y):
                    correct_arrows = {}
                    flipped_arrows = {}
                    for i,arrow in enumerate(transformations[0]):
                        for a in arrow:
                            if any(x != True for x in a):
                                if 'internal electron movement' in a:
                                    n_arrows = len(a)
                                    
                                    for trans in a:
                                        if trans != "internal electron movement":
                                            arrow = trans
                                    correct_arrows[i+1] = [arrow] * n_arrows
                               
                                elif 'arrow_flipped' in a:
                
                                    arrow_flipped = True
                                    flipped_arrows[i+1] = 'arrow flipped'
                                else:
                                    correct_arrows[i+1] = a
                                    break
                    k = list(correct_arrows.keys())
            
                    for x in k:
                        same_arrow = False
                        if all(y == correct_arrows[x][0] for y in correct_arrows[x]):
                            
                            same_arrow = True                 
                        if same_arrow == True:
                            #print(correct_arrows[x])
                            str_tuple = (x,len(correct_arrows[x]),str(correct_arrows[x][0]).replace('[','').replace(']',''))
                        elif same_arrow == False:
                            #print(correct_arrows)
                            str_tuple = (x,len(correct_arrows[x]),str(correct_arrows[x]).replace('[','').replace(']',''))
                        if str_tuple[1] == 1:
                            feedback += "In step %s you should draw %s arrow for a %s. <br>" %str_tuple
                        else:
                            feedback += "In step %s you should draw %s arrows for a %s. <br>" %str_tuple
                f = list(flipped_arrows.keys())        
                if len(f) == 1:
                    feedback += "In step %s make sure the direction the electron flow (curly arrows) is correct. Always start the curved arrow from a lone pair of electrons or a pi-bond.<br>" % list_to_string(f)
                elif len(f) > 1:    
                    feedback += "In steps %s make sure the direction the electron flow (curly arrows) is correct. Always start the curved arrow from a lone pair of electrons or a pi-bond.<br>" % list_to_string(f)
            all_false = all(not any(v) if isinstance(v, list) else not v for v in transformations[1].values())
            if all_false == False:    
                all_steps = transformations[1].keys()
                false_arrow = []
                for x in all_steps:
                    if transformations[1][x] == True:
                        false_arrow.append(x+1)
                if len(false_arrow) == 1:
                    feedback += "In step %s you should reconsider the number of electrons indicated in the electron push (half or full arrow head).<br>" % list_to_string(false_arrow)
                elif len(false_arrow) >1:
                    feedback += "In steps %s you should reconsider the number of electrons indicated in the electron push (half or full arrow head).<br>" % list_to_string(false_arrow)
        if issues["structure"] == False:
            feedback += "<br>Make sure to pay close attention to how you draw the molecules.<br>"
            student_steps = sorted(mol_struc.keys(), key=lambda x: (0 if 'reactant' in x else (2 if 'product' in x else 1)))
            structure_comparison = reaction_steps[0]["molecular_structures"]
            for i,x in enumerate(student_steps):
                if any(v is False for v in mol_struc[x].values()) or (i) in structure_comparison:
                    if x == 'reactants':
                        if len(mol_struc[x].keys()) == 1:
                            feedback += "Check the molecular structure of your reactants for valency errors.<br>"
                        else:
                            feedback += "Check the molecular structure of your reactants for valency errors.<br>"
                    if x == 'products':
                        if len(mol_struc[x].keys()) == 1:
                            feedback += "Check the molecular structure of your products for valency errors.<br>"
                        else:
                            feedback += "Check the molecular structure of your products for valency errors.<br>"
                    if 'intermediate' in x:
                        if sum(1 for value in mol_struc[x].values() if value is False) == 1 and not i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates in step %s for valency errors, as one is not correct.<br>" %str(i + 1)
                        elif sum(1 for value in mol_struc[x].values() if value is False) == 1 and i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s for valency errors. Here, the %s.<br>" %(structure_comparison[i][0][1],str(i + 1),structure_comparison[i][0][0])
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and not i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates in step %s for valency errors, as multiple are not correct.<br>" %str(i + 1)
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates with %s carbons in step %s for valency errors. Here, the %s.<br>" %(structure_comparison[i][0][1],str(i + 1),structure_comparison[i][0][0])
                        elif all(x == True for x in mol_struc[x].values()) and structure_comparison[i]:
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s. Here, the %s.<br>" %(structure_comparison[i][0][1],str(i + 1),structure_comparison[i][0][0])
                            
        if issues["resonance"] == False:
            feedback += "<br>Include a correct resonance structure.<br>"
            incorrect_steps = []
            for key in reaction_steps[5]['resonance']:
                key = key[1].split('_') 
                incorrect_steps.append(int(key[1])+1)
            if len(incorrect_steps) < 2:
                feedback += "In step %s you do not include the correct resonance structure.<br>" % list_to_string(incorrect_steps)
            else:
                feedback += "In steps %s you do not include the correct resonance structure.<br>" % list_to_string(incorrect_steps)

            if reaction_steps[5]['resonance_present'] == False:
                feedback += "Check in the exersice whether you should include a resonance structure."
                
        feedback += "<br>Good luck with the adjustments to your mechanism! Paste the JSON code of your new reaction mechanism in the field above for a new analysis."
        
    else:
        feedback = generate_prompt(openai_key,model,student,reaction_steps,transformations,issues)

    return feedback

def generate_error_mechanism(arrows):
    """ In case the mechanism does not follow the correct overall structure, a feedback is generated here."""

    wrong_arrows = []
    for a in arrows.keys():
        if arrows[a] == False:
            wrong_arrows.append(a)
    
    if len(wrong_arrows) == 1:
        return "The mechanism you submitted does not follow the correct structure.<br><br> You did not draw the %s arrow in the correct direction. Make sure you draw the mechanism horizontally, with reaction arrows separating each step. Resonance structures should be drawn below the main chain, with a vertical two headed resonance arrow. <br><br> Good luck with the adjustments to your mechanism! Paste the JSON code of your new reaction mechanism in the field above for a new analysis." % list_to_string(wrong_arrows)       
    if len(wrong_arrows) > 1:
        return "The mechanism you submitted does not follow the correct structure.<br><br> You did not draw the %s arrows in the correct direction. Make sure you draw the mechanism horizontally, with reaction arrows separating each step. Resonance structures should be drawn below the main chain, with a vertical resonance arrow. <br><br> Good luck with the adjustments to your mechanism! Paste the JSON code of your new reaction mechanism in the field above for a new analysis." % list_to_string(wrong_arrows)       
