
def analysis_feedback(steps,transformations,mol_struc):
    issues = {'global':True,'mechanistic':True,'structure':True,'resonance':True,'rogue':True}

    s_keys = steps['steps_check']['individual_steps'].keys()
    steps_bool = []
    for x in s_keys:
        try:
            steps_bool.append(steps['steps_check']['individual_steps'][x][0][0])
        except:
            steps_bool.append(steps['steps_check']['individual_steps'][x][0])

    if any(False in sublist for sublist in steps_bool):
        if steps['steps_check']["all_steps_match"] == False:
            issues['global'] = False
        elif False in steps['steps_check']["product"].values():
            issues['global'] = False
            
    if steps['secondary_check']['default'] == 'student_default':
        issues['global'] = False
        
    if any(x == False for x in steps['secondary_check']['reaction_arrows'].values()):
        issues['global'] = False
    
    if steps['steps_check']['missing_steps']:
        issues['global'] = False
            
    if transformations and any(x == False for x in transformations[3].values()):
        issues['mechanistic'] = False
        
    if transformations and any(x == True for x in transformations[1].values()):
        issues['mechanistic'] = False

    if any(False in d.values() for d in mol_struc.values()) or len(steps['steps_check']["molecular_structures"].keys()) > 0:
        issues['structure'] = False
        
    if any(False in d.values() for d in mol_struc.values()) and len(steps['steps_check']["molecular_structures"].keys()) > 0:
        issues['structure'] = False
        
    if steps['secondary_check']['resonance'] and any(item[0] != True for item in steps['secondary_check']['resonance']):
        if steps['secondary_check']['resonance'] == ['no_resonance']:
            pass
        else:
            issues['resonance'] = False
    
    if steps['secondary_check']['resonance_present'] == True and steps['secondary_check']['resonance'][0][0] != True:
        issues['resonance'] = False
    
    if steps['rogue_species']:
        issues['rogue'] = False
        
    return issues