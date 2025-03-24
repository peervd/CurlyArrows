
def analysis_feedback(steps,transformations,mol_struc):
    issues = {'global':True,'mechanistic':True,'structure':True}
    
    s_keys = steps[0]['individual_steps'].keys()
    steps_bool = []
    for x in s_keys:
        steps_bool.append(steps[0]['individual_steps'][x][0])

    if any(False in sublist for sublist in steps_bool):
        issues['global'] = False

    if any(x == False or isinstance(x, str) for z in transformations for y in z for x in y):
        issues['mechanistic'] = False
        
    if any(False in d.values() for d in mol_struc.values()):
        issues['structure'] = False

    return issues