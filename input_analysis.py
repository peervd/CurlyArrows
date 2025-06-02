from operate_analysis import analyze

'''

This script can be run to get custom feedback from your openAI account on a 
reaction mechanism. You can draw your reaction mechanism in ChemDoodle. To run
the analysis you have to copy your openai key, the exercise number and the
ChemDoodle drawing as JSON code. 

'''

# Copy the openAI key between quotation marks
# find your key under "specify later"
ai_key = False
 
# integer number which corresponds to the exercise number of the reaction mechanism
exercise =  1  # paste integer number between 1 and 9

# copy json_code from your reaction mechanism in ChemDoodle
json_code = '''
paste json code here

'''


''' now run the script '''
output = analyze(ai_key,exercise,json_code)
print(output)