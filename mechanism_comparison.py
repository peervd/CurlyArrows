from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator,rdFMCS, BRICS
from rdkit.DataStructs import TanimotoSimilarity
from exersices_ORC import load_acid_base
from compare_molecular_structures import compare_molecules
import re

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def sanitize_hypervalent_smiles(smiles):
    """
    Preprocess SMILES strings with hypervalent atoms to make them compatible with RDKit.
    This function implements strategies to handle molecules where atoms exceed normal valence.
    
    Strategies include:
    1. Converting charged carbon centers to neutral forms with appropriate valence
    2. Adding explicit hydrogens where needed
    3. Using wildcard atoms for problematic centers
    
    Parameters:
      smiles (str): Original SMILES string that may contain hypervalent atoms
      
    Returns:
      str: Modified SMILES string that RDKit can parse
    """
    # Strategy 1: Handle hypervalent charged carbons like [C+](C)(C)(C)(C)C
    # Convert patterns like [C+](X)(X)(X)(X)X to [C+](X)(X)(X)X-X
    pattern = r'\[C\+\](\([^)]+\)){5,}'
    while re.search(pattern, smiles):
        match = re.search(pattern, smiles)
        if match:
            start, end = match.span()
            # Extract the excessive bonds and convert to a more reasonable form
            excessive_part = smiles[start:end]
            # Count the number of parentheses groups
            groups = re.findall(r'\([^)]+\)', excessive_part)
            if len(groups) > 4:  # If more than 4 bonds to carbon
                # Keep 4 groups in the charged carbon, move others to separate bonds
                modified = '[C+]' + ''.join(groups[:4])
                remaining = '.'.join([g.strip('()') for g in groups[4:]])
                smiles = smiles[:start] + modified + smiles[end:] + '.' + remaining
    
    # Strategy 2: Handle specific patterns like [H][OH+]C(C)(C)(C)(C)C
    # Convert patterns with too many substituents on carbon
    pattern = r'C(\([^)]+\)){5,}'
    while re.search(pattern, smiles):
        match = re.search(pattern, smiles)
        if match:
            start, end = match.span()
            # Extract the excessive bonds and convert to a more reasonable form
            excessive_part = smiles[start:end]
            # Count the number of parentheses groups
            groups = re.findall(r'\([^)]+\)', excessive_part)
            if len(groups) > 4:  # If more than 4 bonds to carbon
                # Keep carbon with 4 groups, separate others
                modified = 'C' + ''.join(groups[:4])
                remaining = '.'.join([g.strip('()') for g in groups[4:]])
                smiles = smiles[:start] + modified + smiles[end:] + '.' + remaining
    
    # Strategy 3: Replace problematic charged structures with reasonable analogs
    # Example: [C+](C)(C)(C)(C)C → C(C)(C)(C)C
    smiles = re.sub(r'\[C\+\](\([^)]+\)){4,}', r'C\1', smiles)
    
    return smiles

def create_mol_safely(smiles):
    """
    Attempt to create an RDKit molecule from a SMILES string, with fallbacks for problematic structures.
    
    Parameters:
      smiles (str): SMILES string to convert to RDKit molecule
      
    Returns:
      RDMol or None: RDKit molecule object, or None if all attempts fail
    """
    # First try normal parsing
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    
    # If normal parsing fails, try sanitization with various options
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            # Try to sanitize with as many properties as possible
            sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
            return mol
    except:
        pass
    
    # If still failing, try sanitize with more permissive valence checking
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is not None:
            # Try with even more limited sanitization
            sanitize_ops = Chem.SanitizeFlags.SANITIZE_FINDRADICALS | \
                          Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | \
                          Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
            Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
            return mol
    except:
        pass
    
    # Last attempt: preprocess the SMILES to handle hypervalent atoms
    try:
        sanitized_smiles = sanitize_hypervalent_smiles(smiles)
        mol = Chem.MolFromSmiles(sanitized_smiles)
        return mol
    except:
        return None

def compute_substructure_score(mol1, mol2, timeout=1):
    """
    Compute a substructure similarity score based on the maximum common substructure (MCS).
    The score is defined as:
         score = (number of atoms in MCS) / (min(number of atoms in mol1, mol2))
    If no MCS is found, the score is 0.
    """
    res = rdFMCS.FindMCS([mol1, mol2], timeout=timeout)
    if res.canceled or res.numAtoms == 0:
        return 0.0
    return res.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())

def combined_similarity(mol1, mol2, fp1, fp2):
    """
    Combine the Tanimoto similarity of Morgan fingerprints and an MCS-based substructure score.
    Here we take the maximum of the two scores.
    """
    t_sim = TanimotoSimilarity(fp1, fp2)
    s_sim = compute_substructure_score(mol1, mol2)
    return max(t_sim, s_sim)

def get_candidate_fragments(smi, min_atoms=4):
    """
    For a given molecule SMILES, return a set of candidate fragments.
    Always include the parent canonical SMILES.
    If the molecule has at least `min_atoms` heavy atoms,
    add its BRICS fragments (if any).
    """
    candidates = set()
    mol = create_mol_safely(smi)

    if mol is None:
        # Even with sanitization, we couldn't create a valid molecule
        # Add the original SMILES as a last resort
        candidates.add(smi)
        return candidates
    
    try:
        parent = Chem.MolToSmiles(mol, canonical=True)
        candidates.add(parent)
    except:
        # If canonicalization fails, use the original SMILES
        candidates.add(smi)
    
    # Only fragment if the molecule is sufficiently large.
    if mol.GetNumHeavyAtoms() >= min_atoms:
        try:
            frags = BRICS.BRICSDecompose(mol)
            # BRICSDecompose returns a set of fragment SMILES.
            candidates.update(frags)
        except Exception as e:
            pass

    return candidates

def canonicalize_reaction(smiles_reaction, sim_threshold=0.5):
    """
    Convert a reaction SMILES into a canonical form by removing spectator molecules.
    
    This function first fragments each reactant and product (using BRICS) to generate candidate
    substructures. Then, for each parent molecule, if any of its candidate fragments finds a matching
    candidate on the opposite side (with combined similarity >= sim_threshold), the parent is retained.
    
    Molecules that do not find any matching partner are considered spectators and removed.
    
    Parameters:
      smiles_reaction (str): Reaction SMILES with reactants and products separated by '>>'
                             and individual molecules separated by '.'.
      sim_threshold (float): Minimum similarity (0–1) required for a candidate match.
    
    Returns:
      str: Canonical reaction SMILES including only the core transformation.
    """
    try:
        reactants_str, products_str = smiles_reaction.split(">>")
    except Exception as e:
        raise ValueError("Reaction SMILES must contain '>>' separating reactants and products.")

    reactant_smiles = reactants_str.split(".")
    product_smiles = products_str.split(".")
    
    # For each parent molecule, generate candidate fragments.
    reactant_candidates = {}  # key: parent SMILES, value: set of candidate SMILES
    product_candidates = {}
    unsanitized_r_mol = []
    unsanitized_p_mol = []
    for smi in reactant_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            reactant_candidates[smi] = cand
        else:
            unsanitized_r_mol.append(smi)
    for smi in product_smiles:
        cand = get_candidate_fragments(smi)
        if cand:
            product_candidates[smi] = cand
        else:
            unsanitized_p_mol.append(smi)

    # Add the unsanitized molecules back as their own candidates
    # This ensures they are still considered in the reaction
    for smi in unsanitized_r_mol:
        reactant_candidates[smi] = {smi}
    for smi in unsanitized_p_mol:
        product_candidates[smi] = {smi}

    # Create a fingerprint generator.
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    # For each candidate fragment, precompute its RDKit molecule and fingerprint.
    def get_fp(smi):
        mol = create_mol_safely(smi)
        if mol is None:
            return None, None
        try:
            fp = generator.GetFingerprint(mol)
            return mol, fp
        except:
            return None, None
    
    reactant_frag_data = {}  # key: candidate SMILES, value: (mol, fp)
    for parent, cand_set in reactant_candidates.items():
        for cand in cand_set:
            mol, fp = get_fp(cand)
            if mol is not None:
                reactant_frag_data[cand] = (mol, fp)
                
    product_frag_data = {}
    for parent, cand_set in product_candidates.items():
        for cand in cand_set:
            mol, fp = get_fp(cand)
            if mol is not None:
                product_frag_data[cand] = (mol, fp)
    
    # Determine which parent molecules (reactant and product) are "involved".
    involved_reactants = set()
    involved_products = set()
    
    # For each reactant parent, if any of its candidate fragments finds a match in any product candidate, mark it as involved.
    for r_parent, r_cands in reactant_candidates.items():
        for r_cand in r_cands:
            if r_cand not in reactant_frag_data:
                continue
            r_mol, r_fp = reactant_frag_data[r_cand]
            for p_parent, p_cands in product_candidates.items():
                for p_cand in p_cands:
                    if p_cand not in product_frag_data:
                        continue
                    p_mol, p_fp = product_frag_data[p_cand]
                    sim = combined_similarity(r_mol, p_mol, r_fp, p_fp)
                    if sim >= sim_threshold:
                        involved_reactants.add(r_parent)
                        involved_products.add(p_parent)
                        # Once one match is found for this reactant parent, no need to check further.
                        break
                else:
                    continue
                break
    
    # Also check from the product side (for safety)
    for p_parent, p_cands in product_candidates.items():
        for p_cand in p_cands:
            if p_cand not in product_frag_data:
                continue
            p_mol, p_fp = product_frag_data[p_cand]
            for r_parent, r_cands in reactant_candidates.items():
                for r_cand in r_cands:
                    if r_cand not in reactant_frag_data:
                        continue
                    r_mol, r_fp = reactant_frag_data[r_cand]
                    sim = combined_similarity(r_mol, p_mol, r_fp, p_fp)
                    if sim >= sim_threshold:
                        involved_reactants.add(r_parent)
                        involved_products.add(p_parent)
                        break
                else:
                    continue
                break
    
    # String-based fallback for molecules that couldn't be properly converted to RDKit mols
    # Get raw smiles that haven't been classified as involved yet
    reactant_raw_smiles = {smi for smi in reactant_smiles if smi not in involved_reactants}
    product_raw_smiles = {smi for smi in product_smiles if smi not in involved_products}
    
    # This is a simple direct comparison of SMILES, which works for exact matches
    for r_smi in reactant_raw_smiles:
        sanitized_r = sanitize_hypervalent_smiles(r_smi)
        for p_smi in product_raw_smiles:
            sanitized_p = sanitize_hypervalent_smiles(p_smi)
            # Simple string comparison as fallback
            if sanitized_r == sanitized_p:
                involved_reactants.add(r_smi)
                involved_products.add(p_smi)
    
    # Build the canonical reaction SMILES from the parent molecules that are involved.
    core_reactants_smiles = ".".join(sorted(involved_reactants))
    core_products_smiles = ".".join(sorted(involved_products))
    return f"{core_reactants_smiles}>>{core_products_smiles}"


def find_longest_matching_sequences(steps_dict):
    """Find the longest non-overlapping matching subsequences between the model and student's sequences, including shift."""
    steps = steps_dict.keys()
    if all(steps_dict[s][1] == 'not_present' for s in steps):
        non_overlapping_matches = []
        return non_overlapping_matches
    matches = []
    for i,step in enumerate(steps):
        if steps_dict[step][1] == ('not_present'):
            matches.append(('not_present',0))
        else:
            shift = step - steps_dict[step][1]
            if steps_dict[step][0] == [True,True]:
                matches.append((step,shift))
            elif steps_dict[step][0] == [True,False]:
                matches.append((step,shift))
            elif steps_dict[step][0] == [False,True]:
                matches.append((step+1,shift))
                
    # Initialize tracking variables
    current_start = None
    current_end = None
    current_shift = None
    non_overlapping_matches = []

    for i, (step, shift) in enumerate(matches):
        if current_start is None:
            # Start a new sequence
            current_start = step
            current_end = step
            current_shift = shift
        elif step == current_end + 1:
            current_end = step 
        elif step == 'not_present' or step == len(matches)-1:
            if current_start == 'not_present' and step != 'not_present':
                new_sequence = (step, step, shift)
            # Store the completed sequence if it's not already present
            else:
                new_sequence = (current_start, current_end, current_shift)
            non_overlapping_matches.append(new_sequence)
            
            # Start a new sequence
            current_start = None 
            current_end = None
            current_shift = None
    
    if current_start != None and current_end != None and current_start != 'not_present' and current_end != 'not_present':
        only_sequence = (current_start, current_end, current_shift)
        non_overlapping_matches.append(only_sequence)

    return sorted(non_overlapping_matches)

def compare_steps(model_canonical,student_canonical, acid_base):
    individual_comparisons = {}
    molecular_structure = {}
    for st_numb,st_rx in enumerate(student_canonical):
        st_rx = st_rx.split('>>')
        st_r = st_rx[0].split('.')
        st_p = st_rx[1].split('.')
        reaction_found = False

        charge_issues = []
        for m_numb,m_rx in enumerate(model_canonical):
            m_rx = m_rx.split('>>')
            m_r = m_rx[0].split('.')
            m_p = m_rx[1].split('.')
            bool_r = []
            r = []
            bool_p = []
            p = []
            matching_mol = []
            unique_mol = []

            if any(x in m_r for x in st_r):
                for x in st_r:
                    if x in m_r:
                        bool_r.append(True)
                        reactant = True
                        r.append(x)
                        matching_mol.append(x)

                    elif any(x in acid[0] for acid in acid_base) and any(acid[0] in m_r for acid in acid_base):
                        bool_r.append(True)
                        reactant = True
                        r.append(x)
                        matching_mol.append(x)
                    elif any(x in base[1] for base in acid_base) and any(base[1] in m_r for base in acid_base):
                        bool_r.append(True)
                        r.append(x)
                        reactant = True
                        matching_mol.append(x)    
                    else: 
                        bool_r.append(False)
                        unique_mol.append(x)
            if any(x in m_p for x in st_p):
                for x in st_p:
                    if x in m_p:
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
                    elif any(x in acid[0] for acid in acid_base) and any(acid[0] in m_p for acid in acid_base):
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
                    elif any(x in base[1] for base in acid_base) and any(base[1] in m_p for base in acid_base):
                        bool_p.append(True)
                        p.append(x)
                        product = True
                        
                    else: 
                        bool_p.append(False)
            unique_from_model = set(m_r) - set(matching_mol)

            if unique_mol:
                differences = []
                if unique_from_model:
                    for x in unique_mol:
                        comparison = []
                        for y in unique_from_model:
                            comparison.append(compare_molecules(x,y))
                        pri = 5
                        for z in comparison:
                            if z[0].endswith('_1') and pri > 0:
                                diff = [z[0].strip('_1'),z[1]]
                                pri = 0
                            elif z[0].endswith('_2') and pri > 1:
                                diff = [z[0].strip('_2'),z[1]]
                                pri = 1
                            elif z[0].endswith('_3') and pri > 2:
                                diff = [z[0].strip('_3'),z[1]]
                                pri = 2
                            elif z[0].endswith('_4') and pri > 3:
                                diff = [z[0].strip('_4'),z[1]]
                                pri = 3
                            elif z[0].endswith('_0') and pri > 4:
                                diff = [z[0].strip('_0'),z[1]]
                                pri = 4
                            
                        differences.append(diff)
                    molecular_structure[st_numb] = differences
                else:
                    molecular_structure[st_numb] = ['molecule_not_present_in_model',sum(1 for atom in Chem.MolFromSmiles(unique_mol[0]).GetAtoms() if atom.GetSymbol() == 'C')]
    
            if bool_r and bool_p and all(x == True for x in bool_r) and all(x == True for x in bool_p):
                individual_comparisons[st_numb] = ([True,True],m_numb)
                reaction_found = True
                break

            elif bool_r and bool_p and all(x == True for x in bool_r) and any(x == False for x in bool_p):
                max_length = max([len(s) for s in m_r])
                if max_length == max([len(s) for s in r]):
                    individual_comparisons[st_numb] = ([True,False],m_numb)
                    reaction_found = True
                    break

            elif bool_r and bool_p and any(x == False for x in bool_r) and all(x == True for x in bool_p):
                max_length = max([len(s) for s in m_p])
                if max_length == max([len(s) for s in p]):
                    individual_comparisons[st_numb] = ([False,True],m_numb)
                    reaction_found = True
                    break
            

        if reaction_found == False:
            individual_comparisons[st_numb] = ([False,False],'not_present')

    return individual_comparisons,molecular_structure

def count_trues(d):
    return sum(value[0].count(True) for value in d.values())

def compare_reactions(model_list, student_list):
    """
    Compare a model reaction sequence with a student's reaction sequence.
    
    The "individual_steps" key in the output dictionary represents a list of 
    Boolean values. Each entry in this list corresponds to a reaction in the 
    student's list and indicates whether that specific reaction appears in the 
    model list after canonicalization (i.e., after removing spectator 
    molecules)
    
    The "matching_sequences" output contains a list of tuples, where each tuple 
    (i,j) represents a subsequence of the model list that is also present in 
    the student's list. Specifically, it means that the transformations from 
    index i to j (inclusive) in the model list appear somewhere in the student 
    list in the correct order, but possibly shifted in position.

    For example, if "matching_sequences": [(5, 10, 0)], this means that 
    transformations from step 5 to step 10 in the model sequence appear 
    somewhere in the student's list in the correct order. Reactants are always
    indexed as 0.
    """

    student_canonical = [canonicalize_reaction(rx) for rx in student_list]

    acid_base = load_acid_base('list_acid_base.csv')
    index_comp = False
    if isinstance(model_list, list) and all(isinstance(item, list) for item in model_list):
        lst_comparisons = []
        for m_lst in model_list:
            model_canonical = [canonicalize_reaction(rx) for rx in m_lst]
            ind_comp = compare_steps(model_canonical,student_canonical,acid_base)
            lst_comparisons.append(ind_comp)    
        index_comp, individual_comparisons = max(enumerate(lst_comparisons), key=lambda x: count_trues(x[1]))
    
    else:            
        model_canonical = [canonicalize_reaction(rx) for rx in model_list]
        individual_comparisons = compare_steps(model_canonical, student_canonical,acid_base)

    matching_subsequences = find_longest_matching_sequences(individual_comparisons[0])

    return {"individual_steps": individual_comparisons[0],"matching_sequences": matching_subsequences,
            "molecular_structures":individual_comparisons[1],"index_resonance": index_comp}