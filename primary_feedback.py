
def analysis_feedback(steps,transformations,mol_struc):
    issues = {'global':True,'mechanistic':True,'structure':True,'resonance':True}
    
    s_keys = steps[0]['individual_steps'].keys()
    steps_bool = []
    for x in s_keys:
        try:
            steps_bool.append(steps[0]['individual_steps'][x][0][0])
        except:
            steps_bool.append(steps[0]['individual_steps'][x][0])
    
    if any(False in sublist for sublist in steps_bool):
        issues['global'] = False
    
    if steps[-1]['default'] == 'student_default':
        issues['global'] = False
        
    if any(x == False for x in steps[-1]['reaction_arrows'].values()):
        issues['global'] = False
        
    if transformations and any(x == False or isinstance(x, str) for z in transformations[0] for y in z for x in y):
        issues['mechanistic'] = False
        
    if transformations and any(x == True for x in transformations[1].values()):
        issues['mechanistic'] = False
    
    if any(False in d.values() for d in mol_struc.values()) or len(steps[0]["molecular_structures"].keys()) > 0:
        issues['structure'] = False
        
    if any(False in d.values() for d in mol_struc.values()) and len(steps[0]["molecular_structures"].keys()) > 0:
        issues['structure'] = False
        
    if steps[5]['resonance'] and any(item[0] != True for item in steps[5]['resonance']):
        if steps[5]['resonance'] == ['no_resonance']:
            pass
        else:
            issues['resonance'] = False
    
    if steps[5]['resonance_present'] == True and steps[5]['resonance'] != ["no_resonance"]:
        issues['resonance'] = False
        
    return issues