def generate_prompt(openai_key, model_reaction, student, reaction_steps, transformations, issues):
    return "Function to be written, based on choise for AI model."

def generate_feedback(openai_key,model,student,reaction_steps,transformations,issues):
    """ Function to generate feedback based on the analysis of the student reaction mechanism and the selected exersice. """
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
            for i,x in enumerate(reaction_steps[0]['individual_steps']):
                if x == False:
                    incorrect_steps.append(i+1)
            feedback += "Your reaction mechanism does not follow the correct steps:<br>"
            if len(incorrect_steps) == 1:
                feedback += "Step %s does not correspond with the model answer.<br><br>" % str(incorrect_steps).replace('[','').replace(']','')
            else:
                feedback += "Steps %s do not correspond with the model answer.<br><br>" % str(incorrect_steps).replace('[','').replace(']','')
    
        if issues["mechanistic"] == False: 
            correct_arrows = {}
            for i,x in enumerate(transformations):
                if any(x != True for x in x[0]):
                    correct_arrows[i+1] = x[0]
            feedback += "Please pay attention to the electron flow (curly arrows):<br>"
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
                    feedback += "In step %s you should draw %s arrow for a %s. <br><br>" %str_tuple
                else:
                    feedback += "In step %s you should draw %s arrows for a %s. <br><br>" %str_tuple
            feedback += "Good luck with the adjustments to your mechanism! Paste the json code of your new reaction mechanism in the field above for a new analysis."
        if issues["structure"] == False:
            feedback += "Check the molecular structure for valency and charge errors.<br><br>"

    else:
        feedback = generate_prompt(openai_key,model,student,reaction_steps,transformations,issues)

    return feedback
