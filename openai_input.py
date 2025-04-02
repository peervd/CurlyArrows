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

def list_to_string(lst):
    if len(lst) > 1:
        return ', '.join(map(str, lst[:-1])) + ' and ' + str(lst[-1])
    return ''.join(map(str, lst))

def generate_feedback(openai_key,model,student,reaction_steps,transformations,mol_struc,issues):
    if reaction_steps[3] == False:
        return "Make sure you select the correct exercise. You do not start your mechanism with the correct reactants"

    """Generates structured feedback based on analysis."""
    if openai_key == False:
        feedback = "Your input was carefully analyzed and here is your feedback without openAI analysis.<br><br>"
        issue_list = list(issues.values())
        if all(x == True for x in issue_list):
            return "Your reaction mechanism is correct."
        
        if issues["global"] == False:
            incorrect_steps = []
            s_keys = reaction_steps[0]['individual_steps'].keys()
            steps_bool = []
            for x in s_keys:
                steps_bool.append(reaction_steps[0]['individual_steps'][x][0])
            for i,x in enumerate(steps_bool):
                if x[0] == False and i == 0:
                    incorrect_steps.append(i + 1)
                if x[1] == False:
                    if (i + 2) not in incorrect_steps:
                        incorrect_steps.append(i + 2)
                
            feedback += "The reaction mechanism does not follow the correct steps:<br>"
            if len(incorrect_steps) == 1:
                feedback += "Step %s does not correspond with the model answer.<br>" % list_to_string(incorrect_steps)
            else:
                feedback += "Steps %s do not correspond with the model answer.<br>" % list_to_string(incorrect_steps)
    
        if issues["mechanistic"] == False: 
            correct_arrows = {}
            flipped_arrows = {}
            for i,x in enumerate(transformations):
                if any(x != True for x in x[0]):
                    
                    if 'internal electron movement' in x[0]:
                        n_arrows = len(x[0])
                        for trans in x[0]:
                            if trans != "internal electron movement":
                                arrow = trans
                        correct_arrows[i+1] = [arrow] * n_arrows
                    
                    elif 'arrow_flipped' in x[1]:
                        arrow_flipped = True
                        flipped_arrows[i+1] = 'arrow flipped'
                    else:
                        correct_arrows[i+1] = x[0]

            feedback += "<br>Please pay attention to the electron flow (curly arrows):<br>"
            k = list(correct_arrows.keys())

            for x in k:
                same_arrow = False
                if all(y == correct_arrows[x][0] for y in correct_arrows[x]):
                    
                    same_arrow = True                 
                if same_arrow == True:
                    str_tuple = (x,len(correct_arrows[x]),str(correct_arrows[x][0]).replace('[','').replace(']',''))
                elif same_arrow == False:
                    str_tuple = (x,len(correct_arrows[x]),str(correct_arrows[x]).replace('[','').replace(']',''))
                if str_tuple[1] == 1:
                    feedback += "In step %s you should draw %s arrow for a %s. <br>" %str_tuple
                else:
                    feedback += "In step %s you should draw %s arrows for a %s. <br>" %str_tuple
            f = list(flipped_arrows.keys())        
            if len(f) == 1:
                feedback += "In step %s make sure the direction the electron flow is correct. Always start the curved arrow from a lone pair of electrons or a pi-bond." % list_to_string(f)
            elif len(f) > 1:    
                feedback += "In steps %s make sure the direction the electron flow is correct. Always start the curved arrow from a lone pair of electrons or a pi-bond." % list_to_string(f)
            
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
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s for valency errors. Here, the %s.<br>" %(structure_comparison[i][0][1],int(re.search(r'\d+$', x).group()) +1,structure_comparison[i][0][0])
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and not i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates in step %s for valency errors, as multiple are not correct.<br>" %int(re.search(r'\d+$', x).group())+1 
                        elif sum(1 for value in mol_struc[x].values() if value is False) > 1 and i in structure_comparison:
                            feedback += "Check the molecular structure of your intermediates with %s carbons in step %s for valency errors. Here, the %s.<br>" %(structure_comparison[i][0][1],int(re.search(r'\d+$', x).group()) +1,structure_comparison[i][0][0])
                        elif all(x == True for x in mol_struc[x].values()) and structure_comparison[i]:
                            feedback += "Check the molecular structure of your intermediate with %s carbons in step %s. Here, the %s.<br>" %(structure_comparison[i][0][1],int(re.search(r'\d+$', x).group()) +1,structure_comparison[i][0][0])
                            
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
                
        feedback += "<br>Good luck with the adjustments to your mechanism! Paste the json code of your new reaction mechanism in the field above for a new analysis."
        
    else:
        feedback = generate_prompt(openai_key,model,student,reaction_steps,transformations,issues)

    return feedback