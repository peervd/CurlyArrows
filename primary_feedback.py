
def analysis_feedback(steps,transformations,mol_struc):
    issues = {'global':True,'mechanistic':True,'structure':True}

    if any(x is False for x in steps[0]['individual_steps']):
        issues['global'] = False

    if any(x == False or isinstance(x, str) for z in transformations for y in z for x in y):
        issues['mechanistic'] = False
        
    if any(False in d.values() for d in mol_struc.values()):
        issues['structure'] = False

    return issues