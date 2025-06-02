from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdFMCS, BRICS
from rdkit.DataStructs import TanimotoSimilarity
import re
RDLogger.DisableLog('rdApp.*')

def getSubmolRadAtom(mol, index, radius):
    """
    Extract a submol around a specific atom with a given radius and
    return the SMILES and the index of the original atom in the submol.
    """
    #print(f"Processing atom {index} with radius {radius}")
    
    # Instead of using FindAtomEnvironmentOfRadiusN, let's manually build
    # the environment by breadth-first traversal from the central atom
    
    # Initialize sets to track atoms and bonds
    atoms_to_include = set([index])  # Start with the central atom
    bonds_to_include = set()
    
    # BFS traversal
    current_radius = 0
    frontier = [index]
    visited = set([index])
    
    while current_radius < radius:
        next_frontier = []
        for atom_idx in frontier:
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                # Get the other atom in this bond
                other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
                
                # Add the bond
                bonds_to_include.add(bond.GetIdx())
                
                # Add the other atom
                atoms_to_include.add(other_atom_idx)
                
                # Add to next frontier if not visited
                if other_atom_idx not in visited:
                    visited.add(other_atom_idx)
                    next_frontier.append(other_atom_idx)
        
        frontier = next_frontier
        current_radius += 1
    
    #print(f"Manually found atoms: {atoms_to_include}")
    #print(f"Manually found bonds: {bonds_to_include}")
    
    # Map old atom indices to new ones
    atom_map = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted(atoms_to_include))}
    
    # Convert bonds_to_include from a set to a list for PathToSubmol
    bonds_list = list(bonds_to_include)
    
    # Create the submol using the list of bonds
    submol = Chem.PathToSubmol(mol, bonds_list, atomMap=atom_map)
    
    # Convert to SMILES
    smiles = Chem.MolToSmiles(submol)
    # Get new index of original atom
    new_index = atom_map[index]
    return smiles, new_index

def getSubmolRadBond(mol, index, radius):
    """
    Extract a submol around a specific bond with a given radius
    
    Parameters:
    mol (Chem.Mol): RDKit molecule
    index (tuple): Tuple of atom indices defining the bond
    radius (int): Radius of bond environment to extract
    
    Returns:
    tuple: (SMILES representation of the extracted submol, 
           tuple with 1-based atom indices in the SMILES string for the bond atoms)
    """

    
    # Ensure index is a tuple of two atom indices
    if not isinstance(index, tuple) or len(index) != 2:
        raise ValueError("Index must be a tuple of two atom indices")
    
    # Find atom environments for both atoms in the bond
    envs = []
    for ind in index:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, ind)
        envs.append(env)
    
    # Combine environments
    combined_env = set()
    for env in envs:
        combined_env.update(env)
    
    # Prepare atom map for submol extraction
    amap = {}
    
    # Extract submol using combined environment
    submol = Chem.PathToSubmol(mol, list(combined_env), atomMap=amap)
    
    # Get the new indices of the two bond atoms in the submolecule
    new_idx1 = amap.get(index[0])
    new_idx2 = amap.get(index[1])
    
    # Get clean SMILES without atom mapping - this is what we'll return
    submol_smiles = Chem.MolToSmiles(submol)
    
    # Now we need to figure out the atom order in the SMILES string
    # We'll use a different approach: generate SMILES with atom mapping and parse it manually
    
    # Create a copy for mapping
    submol_copy = Chem.Mol(submol)
    
    # Set atom mapping numbers
    for atom in submol_copy.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)  # 1-based mapping
    
    # Generate mapped SMILES
    mapped_smiles = Chem.MolToSmiles(submol_copy)
    
    # Parse the mapped SMILES manually to extract atom order
    # Look for patterns like :1], :2], etc. and track their order of appearance
    atom_order = []
    i = 0
    while i < len(mapped_smiles):
        if mapped_smiles[i] == ':':
            # Found an atom mapping, extract the number
            j = i + 1
            while j < len(mapped_smiles) and mapped_smiles[j].isdigit():
                j += 1
            if j < len(mapped_smiles) and mapped_smiles[j] == ']':
                atom_map_num = int(mapped_smiles[i+1:j])
                atom_order.append(atom_map_num - 1)  # Convert back to 0-based
        i += 1
    
    # Create mapping from submol atom index to SMILES position (1-based)
    smiles_position = {}
    for smiles_idx, submol_idx in enumerate(atom_order):
        smiles_position[submol_idx] = smiles_idx + 1  # 1-based position in SMILES
    
    # Get the positions for our bond atoms
    pos1 = smiles_position.get(new_idx1, -1) if new_idx1 is not None else -1
    pos2 = smiles_position.get(new_idx2, -1) if new_idx2 is not None else -1
    
    return submol_smiles, (pos1, pos2)

def getRadBond(mol, index):
    """
    Extract SMILES for a specific bond
    
    Parameters:
    mol (Chem.Mol): RDKit molecule
    index (tuple): Tuple of two atom indices defining the bond
    
    Returns:
    str: SMILES representation of the bond
    """
    # Ensure index is a tuple of two atom indices
    if not isinstance(index, tuple) or len(index) != 2:
        raise ValueError("Index must be a tuple of two atom indices")
    
    # Extract bond SMILES
    submol = Chem.MolFragmentToSmiles(mol, list(index))
    return submol

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


def compare_molecules(smiles1, smiles2):
    """
    Compare two molecules represented by SMILES strings to determine their relationship.
    
    Parameters:
    -----------
    smiles1 : str
        SMILES string of the first molecule
    smiles2 : str
        SMILES string of the second molecule
        
    Returns:
    --------
    str
        A string describing the relationship between the molecules:
        - "Same structure but different charges"
        - "First molecule is a substructure of the second molecule"
        - "Second molecule is a substructure of the first molecule"
        - "Molecules are completely different"
        - "Invalid SMILES string(s)" if parsing fails
    """
    # Try to create RDKit molecule objects from SMILES
    try:
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 is None or mol2 is None:
            return "Invalid SMILES string(s)_0"
        
        # Count carbon atoms in the first molecule
        carbon_count = sum(1 for atom in mol1.GetAtoms() if atom.GetSymbol() == 'C')
        
        # Get formal charges for each atom in both molecules
        charges1 = [atom.GetFormalCharge() for atom in mol1.GetAtoms()]
        charges2 = [atom.GetFormalCharge() for atom in mol2.GetAtoms()]
        
        # Remove all charges to check structural identity
        mol1_neutral = Chem.Mol(mol1)
        mol2_neutral = Chem.Mol(mol2)
        
        for atom in mol1_neutral.GetAtoms():
            atom.SetFormalCharge(0)
        
        for atom in mol2_neutral.GetAtoms():
            atom.SetFormalCharge(0)
            
        # Check if the neutral structures are identical
        same_structure = (Chem.MolToSmiles(mol1_neutral) == Chem.MolToSmiles(mol2_neutral))
        
        if same_structure and charges1 != charges2:
            return ["molecule does not have the correct charge_1",carbon_count]
        
        # Check substructure relationships
        if mol1.HasSubstructMatch(mol2):
            if mol2.HasSubstructMatch(mol1):
                # If they match each other but aren't the same structure with different charges,
                # then they might be equivalent representations or have different stereo/isotope info
                return ["molecule should have different stereochemistry/isotopes_2",carbon_count]
            else:
                return ["molecule contains an incorrectly added structure_3",carbon_count]
        elif mol2.HasSubstructMatch(mol1):
            return ["molecule misses part of its structure_4",carbon_count]
        else:
            return ["molecule does not match the model answer_0",carbon_count]
            
    except Exception as e:
        return f"Error comparing molecules: {str(e)}"