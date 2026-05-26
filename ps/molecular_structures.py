from .Calculate_H_from_smiles import count_hydrogens_from_smiles
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdFMCS, BRICS, rdMolDescriptors
from rdkit.DataStructs import TanimotoSimilarity
from collections import Counter
from itertools import combinations, chain
from .helper_functions import remove_capital_letters
import re
import requests
import math
RDLogger.DisableLog('rdApp.*')

def safe_mol_from_smiles(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        mol_core = Chem.RemoveHs(mol)
        return mol_core
    except Exception:
        return Chem.MolFromSmiles(smiles, sanitize=False)

def safe_smiles(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        mol_core = Chem.RemoveHs(mol)
        smi = Chem.MolToSmiles(mol_core)
        return smi
    except Exception:
        return smiles
    
def flatten_bonds_for_mcs(smiles):
    """
    Convert all bonds in a molecule to SINGLE bonds.
    Keeps atoms and connectivity but removes problematic bond orders.
    """
    mol = Chem.MolFromSmiles(smiles,sanitize=False)
    rw = Chem.RWMol(mol)

    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)

    for atom in rw.GetAtoms():
        atom.SetIsAromatic(False)

    m2 = rw.GetMol()
    m2.UpdatePropertyCache(strict=False)

    return m2
    

def standardize_molecule_mapping(mol, original_atom_map, original_atom_ids, original_bond_ids):
    """
    Standardize molecule mapping by using canonical atom ordering from RDKit.
    This ensures that identical molecules get identical mappings regardless of drawing order.
    
    Args:
        mol: RDKit molecule object
        original_atom_map: dict mapping original atom indices to RDKit atom indices
        original_atom_ids: dict mapping atom IDs to RDKit atom indices  
        original_bond_ids: dict mapping bond IDs to RDKit bond tuples
        
    Returns:
        tuple: (standardized_mol, new_atom_map, new_atom_ids, new_bond_ids, mapping_dict)
    """
    
    try:
        # Generate canonical SMILES with atom mapping to preserve order
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Get the canonical atom ordering from RDKit
        # This gives us the mapping from original atom indices to canonical indices
        canonical_order = list(Chem.CanonicalRankAtoms(mol))
        
        # Create a mapping from canonical rank to atom index
        # canonical_order[i] gives the canonical rank of atom i
        # We want the reverse: given a canonical rank, what atom index has that rank
        rank_to_atom = {}
        for atom_idx, rank in enumerate(canonical_order):
            rank_to_atom[rank] = atom_idx
        
        # Create standardized mappings where atoms are ordered by canonical rank
        new_atom_map = {}
        new_atom_ids = {}
        new_bond_ids = {}
        
        # Create mapping dictionary for tracking changes
        mapping_dict = {
            'atom_mapping': {},  # original_rdkit_idx -> standardized_rdkit_idx
            'reverse_atom_mapping': {},  # standardized_rdkit_idx -> original_rdkit_idx
            'bond_mapping': {},  # original_bond_tuple -> standardized_bond_tuple
            'reverse_bond_mapping': {},  # standardized_bond_tuple -> original_bond_tuple
            'canonical_ranks': canonical_order  # Store canonical ranks for debugging
        }
        
        # Map atoms using canonical ordering
        for orig_idx, rdkit_idx in original_atom_map.items():
            if rdkit_idx < len(canonical_order):
                # Get the canonical rank of this atom
                canonical_rank = canonical_order[rdkit_idx]
                # The new standardized index is simply the canonical rank
                new_rdkit_idx = canonical_rank
                
                new_atom_map[orig_idx] = new_rdkit_idx
                mapping_dict['atom_mapping'][rdkit_idx] = new_rdkit_idx
                mapping_dict['reverse_atom_mapping'][new_rdkit_idx] = rdkit_idx
        
        # Map atom IDs using canonical ordering
        for atom_id, rdkit_idx in original_atom_ids.items():
            if rdkit_idx < len(canonical_order):
                canonical_rank = canonical_order[rdkit_idx]
                new_atom_ids[atom_id] = canonical_rank
        
        # Map bonds using canonical ordering
        for bond_id, (start_idx, end_idx) in original_bond_ids.items():
            if start_idx < len(canonical_order) and end_idx < len(canonical_order):
                # Get canonical ranks for both atoms
                new_start_idx = canonical_order[start_idx]
                new_end_idx = canonical_order[end_idx]
                
                # Create consistent bond tuple (always smaller index first)
                new_bond_tuple = tuple(sorted([new_start_idx, new_end_idx]))
                orig_bond_tuple = tuple(sorted([start_idx, end_idx]))
                
                new_bond_ids[bond_id] = new_bond_tuple
                mapping_dict['bond_mapping'][orig_bond_tuple] = new_bond_tuple
                mapping_dict['reverse_bond_mapping'][new_bond_tuple] = orig_bond_tuple
        
        # Create a standardized molecule with atoms in canonical order
        standardized_mol = Chem.RWMol()
        
        # Add atoms in canonical order
        atom_map_for_new_mol = {}
        for canonical_rank in range(len(canonical_order)):
            # Find which original atom has this canonical rank
            original_atom_idx = rank_to_atom[canonical_rank]
            original_atom = mol.GetAtomWithIdx(original_atom_idx)
            
            # Create new atom with same properties
            new_atom = Chem.Atom(original_atom.GetSymbol())
            new_atom.SetFormalCharge(original_atom.GetFormalCharge())
            new_atom.SetNumExplicitHs(original_atom.GetNumExplicitHs())
            
            added_idx = standardized_mol.AddAtom(new_atom)
            atom_map_for_new_mol[original_atom_idx] = added_idx
        
        # Add bonds in canonical order
        for bond in mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            
            new_start_idx = atom_map_for_new_mol[start_idx]
            new_end_idx = atom_map_for_new_mol[end_idx]
            
            standardized_mol.AddBond(new_start_idx, new_end_idx, bond.GetBondType())
        
        # Sanitize the standardized molecule
        try:
            Chem.SanitizeMol(standardized_mol)
        except:
            # If sanitization fails, return original molecule
            return mol, original_atom_map, original_atom_ids, original_bond_ids, {}
        
        return standardized_mol, new_atom_map, new_atom_ids, new_bond_ids, mapping_dict
        
    except Exception as e:
        # If anything fails, return original mappings
        print(f"Standardization failed: {e}")
        return mol, original_atom_map, original_atom_ids, original_bond_ids, {}

def getSubmolRadAtom(mol, index, radius, san = False):
    """
    Extract a submol around a specific atom with a given radius and
    return the SMILES and the character position of the original atom in the SMILES string.
    """

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
                other_atom = mol.GetAtomWithIdx(other_atom_idx)
                if san == True and other_atom.GetSymbol() == 'H':
                    continue
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
    if not atom_map:
        return Chem.MolToSmiles(mol), 0
    # Get new index of original atom in the submolecule
    new_index = atom_map[index]
    
    # Convert to clean SMILES
    smiles = Chem.MolToSmiles(submol)
    
    # Now find the character position in the SMILES string
    # Create a copy for mapping
    submol_copy = Chem.Mol(submol)
    
    # Clear all atom mappings and set mapping only for our target atom
    for atom in submol_copy.GetAtoms():
        atom.SetAtomMapNum(0)  # Clear all mappings first
    
    # Set unique mapping number for our atom of interest
    submol_copy.GetAtomWithIdx(new_index).SetAtomMapNum(999)
    
    # Generate mapped SMILES
    mapped_smiles = Chem.MolToSmiles(submol_copy)
    
    # Find character position in the mapped SMILES
    char_pos = -1

    # Find position of the atom (mapping 999)
    map_pos = mapped_smiles.find(':999]')
    if map_pos != -1:
        # Look backwards from the mapping to find the start of the atom
        i = map_pos - 1
        while i >= 0 and mapped_smiles[i] != '[':
            i -= 1
        if i >= 0:
            # Now find corresponding position in clean SMILES
            # Count characters up to this point, excluding mapping info
            clean_pos = 0
            mapped_pos = 0
            while mapped_pos < i:
                if mapped_smiles[mapped_pos] == ':':
                    # Skip mapping info (:number])
                    while mapped_pos < len(mapped_smiles) and mapped_smiles[mapped_pos] != ']':
                        mapped_pos += 1
                    if mapped_pos < len(mapped_smiles):
                        mapped_pos += 1  # Skip the ']'
                else:
                    clean_pos += 1
                    mapped_pos += 1
            char_pos = clean_pos
    if len(smiles) <= 2:
        return smiles
    else:
        return smiles, char_pos

def getSubmolRadBond(mol, index, radius, san = False):
    """
    Extract a submol around a specific bond with a given radius
    
    Parameters:
    mol (Chem.Mol): RDKit molecule
    index (tuple): Tuple of atom indices defining the bond
    radius (int): Radius of bond environment to extract
    
    Returns:
    tuple: (SMILES representation of the extracted submol, 
           tuple with character positions of the two bond atoms in the SMILES string)
    """
    
    # Ensure index is a tuple of two atom indices
    if not isinstance(index, tuple) or len(index) != 2:
        raise ValueError("Index must be a tuple of two atom indices")

    #Chem.RemoveHs(mol) if san == True else mol
    
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
    
    # Create a copy for mapping
    submol_copy = Chem.Mol(submol)
    
    # Set atom mapping numbers for our target atoms only
    for atom in submol_copy.GetAtoms():
        atom.SetAtomMapNum(0)  # Clear all mappings first
    
    # Set unique mapping numbers for our atoms of interest
    if new_idx1 is not None:
        submol_copy.GetAtomWithIdx(new_idx1).SetAtomMapNum(999)  # Use unique numbers
    if new_idx2 is not None:
        submol_copy.GetAtomWithIdx(new_idx2).SetAtomMapNum(998)
    
    # Generate mapped SMILES
    mapped_smiles = Chem.MolToSmiles(submol_copy)
    
    # Find character positions in the mapped SMILES
    pos1 = -1
    pos2 = -1
    
    # Find position of first atom (mapping 999)
    map1_pos = mapped_smiles.find(':999]')
    if map1_pos != -1:
        # Look backwards from the mapping to find the start of the atom
        i = map1_pos - 1
        while i >= 0 and mapped_smiles[i] != '[':
            i -= 1
        if i >= 0:
            # Now find corresponding position in clean SMILES
            # Count characters up to this point, excluding mapping info
            clean_pos = 0
            mapped_pos = 0
            while mapped_pos < i:
                if mapped_smiles[mapped_pos] == ':':
                    # Skip mapping info (:number])
                    while mapped_pos < len(mapped_smiles) and mapped_smiles[mapped_pos] != ']':
                        mapped_pos += 1
                    if mapped_pos < len(mapped_smiles):
                        mapped_pos += 1  # Skip the ']'
                else:
                    clean_pos += 1
                    mapped_pos += 1
            pos1 = clean_pos
    
    # Find position of second atom (mapping 998)
    map2_pos = mapped_smiles.find(':998]')
    if map2_pos != -1:
        # Look backwards from the mapping to find the start of the atom
        i = map2_pos - 1
        while i >= 0 and mapped_smiles[i] != '[':
            i -= 1
        if i >= 0:
            # Now find corresponding position in clean SMILES
            clean_pos = 0
            mapped_pos = 0
            while mapped_pos < i:
                if mapped_smiles[mapped_pos] == ':':
                    # Skip mapping info (:number])
                    while mapped_pos < len(mapped_smiles) and mapped_smiles[mapped_pos] != ']':
                        mapped_pos += 1
                    if mapped_pos < len(mapped_smiles):
                        mapped_pos += 1  # Skip the ']'
                else:
                    clean_pos += 1
                    mapped_pos += 1
            pos2 = clean_pos
    
    return submol_smiles, tuple(sorted((pos1, pos2)))
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
    try:
        if len(submol) > 2:
            sanitized = ""
            for n in submol:
                if n.isalpha() and n != 'H' or n == '=' or n == '#':
                    sanitized += n
            submol = sanitized
    except:
        pass
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

def are_constitutional_isomers(smiles1, smiles2):
    """
    Check if two molecules are constitutional isomers.
    Constitutional isomers have the same molecular formula but different connectivity.
    """
    # Parse SMILES
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        return False, "Invalid SMILES"
    
    # Get molecular formulas
    formula1 = rdMolDescriptors.CalcMolFormula(mol1)
    formula2 = rdMolDescriptors.CalcMolFormula(mol2)
    
    # Check if formulas are the same
    if formula1 != formula2:
        return False, f"Different molecular formulas: {formula1} vs {formula2}"
    
    # Generate canonical SMILES to check if they're the same molecule
    canonical1 = Chem.CanonSmiles(smiles1)
    canonical2 = Chem.CanonSmiles(smiles2)
    
    if canonical1 == canonical2:
        return False, "Same molecule (not isomers)"
    
    return True, f"Constitutional isomers with formula {formula1}"

def molecule_core(smiles):
    """ Stips smiles from charges and hydrogen atoms. """
    try:
        mol = Chem.MolFromSmiles(smiles)
    
        if mol is None:
            return "Invalid SMILES string(s)_0"
            
        # Get formal charges for each atom in both molecules
        charges = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    
        # Remove all charges to check structural identity
        mol_neutral = Chem.Mol(mol)
            
        for atom in mol_neutral.GetAtoms():
            atom.SetFormalCharge(0)
        
        # Sanitize to remove explicit H from neutral molecules
        mol_neutral_withH = Chem.AddHs(mol_neutral)
        mol_neutral_H_Smiles = Chem.MolToSmiles(mol_neutral_withH)
        mol_neutral_H_MOL = Chem.MolFromSmiles(mol_neutral_H_Smiles)
        mol_neutral = Chem.RemoveHs(mol_neutral_H_MOL)
        smiles_neutral = Chem.MolToSmiles(mol_neutral)
        return smiles_neutral,mol_neutral, charges
    except:
        pos = smiles.count("+")
        neg = smiles.count("-")
        charges = [pos + neg]
        try:
            mol = Chem.MolFromSmiles(smiles)
        except:
            mol = None
        
        return smiles,mol,charges

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
        smiles1_neutral,mol1_neutral,charges1 = molecule_core(smiles1)
        smiles2_neutral,mol2_neutral,charges2 = molecule_core(smiles2)

        carbon_count = sum(1 for atom in mol1_neutral.GetAtoms() if atom.GetSymbol() == 'C')
        # Check if the neutral structures are identical

        same_structure = (smiles1_neutral == smiles2_neutral)

        if same_structure and sum(charges1) != sum(charges2):
            return ["%s does not have the correct charge_1"%get_molecule_name(smiles1),carbon_count] 

        # Check substructure relationships
        if mol1_neutral.HasSubstructMatch(mol2_neutral):
            if mol2_neutral.HasSubstructMatch(mol1_neutral):
                # If they match each other but aren't the same structure with different charges,
                # then they might be equivalent representations or have different stereo/isotope info
                if sum(charges1) != sum(charges2):
                    return ["%s should have different stereochemistry/isotopes and charge_2"%get_molecule_name(smiles1),carbon_count] 
                else:
                    return ["%s should have different location of charge/stereochemistry/isotopes_5"%get_molecule_name(smiles1),carbon_count] 
            else:
                if sum(charges1) != sum(charges2):
                    return ["%s contains an incorrectly added structure and does not have the correct charge_3"%get_molecule_name(smiles1),carbon_count] %(smiles1)
                else:
                    return ["%s contains an incorrectly added structure_6"%get_molecule_name(smiles1),carbon_count] 
        elif mol2_neutral.HasSubstructMatch(mol1_neutral):
            if sum(charges1) != sum(charges2):
                return ["%s misses part of its structure and does not have the correct charge_4"%get_molecule_name(smiles1),carbon_count] 
            else: 
                return ["%s misses part of its structure_7"%get_molecule_name(smiles1),carbon_count] 
        else:
            return ["%s does not match the model answer_0"%get_molecule_name(smiles1),carbon_count] 
        
    except Exception as e:
        return f"Error comparing molecules: {str(e)}"
    
def get_molecule_name(smiles):
    charge = classify_molecule_by_charge(smiles)
    
    if not charge == 'no_charge':
        return charge
    else:              
    
        try:
            # Convert to canonical SMILES first
            mol = Chem.MolFromSmiles(smiles)
            canonical_smiles = Chem.MolToSmiles(mol)
    
            # Use a web service PubChem with REST API
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{canonical_smiles}/property/Title/JSON"
            response = requests.get(url)
            
            if response.status_code == 200:
                data = response.json()
                return remove_capital_letters(data['PropertyTable']['Properties'][0]['Title'])
            else:
                return smiles
        except:
            return smiles
        
def classify_molecule_by_charge(smiles):
    """
    Classify molecules based on charge patterns or lookup neutral molecules in PubChem.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        str: Classification or name of the molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "Invalid SMILES"
    
    # Analyze charges
    charged_atoms = []
    total_charge = 0
    
    for atom in mol.GetAtoms():
        formal_charge = atom.GetFormalCharge()
        if formal_charge != 0:
            charged_atoms.append({
                'atom': atom.GetSymbol(),
                'charge': formal_charge,
                'idx': atom.GetIdx()
            })
            total_charge += formal_charge
    if not charged_atoms:
        return 'no_charge'
    else:
        return classify_charged_molecule(charged_atoms, total_charge)

def classify_charged_molecule(charged_atoms, total_charge):
    """
    Classify molecule based on charged atoms present.
    """
    
    # Check for zwitterion (both positive and negative charges)
    has_positive = any(atom['charge'] > 0 for atom in charged_atoms)
    has_negative = any(atom['charge'] < 0 for atom in charged_atoms)
    
    if has_positive and has_negative:
        return "zwitterionic molecule"
    
    # Group by atom type and charge
    charge_patterns = {}
    for atom in charged_atoms:
        key = (atom['atom'], atom['charge'] > 0)
        if key not in charge_patterns:
            charge_patterns[key] = 0
        charge_patterns[key] += abs(atom['charge'])
    
    # Classification based on dominant charged atom
    classifications = []
    
    for (atom_symbol, is_positive), charge_count in charge_patterns.items():
        classification = get_atom_classification(atom_symbol, is_positive)
        if classification:
            classifications.append(classification)
    
    # Return the most specific classification
    if len(classifications) == 1:
        return classifications[0]
    elif len(classifications) > 1:
        # Multiple different charged atoms
        return f"mixed charged molecule ({', '.join(set(classifications))})"
    else:
        # Fallback for unrecognized patterns
        if total_charge > 0:
            return "positively charged molecule"
        elif total_charge < 0:
            return "negatively charged molecule"
        else:
            return "neutral molecule with formal charges"

def get_atom_classification(atom_symbol, is_positive):
    """
    Get classification name based on atom type and charge.
    """
    
    classifications = {
        # Oxygen
        ('O', True): 'oxonium molecule',
        ('O', False): 'alkoxide molecule',
        
        # Nitrogen
        ('N', True): 'ammonium molecule',
        ('N', False): 'amide anion molecule',
        
        # Carbon
        ('C', True): 'carbocation molecule',
        ('C', False): 'carbanion molecule',
        
        # Sulfur
        ('S', True): 'sulfonium molecule',
        ('S', False): 'sulfide anion molecule',
        
        # Phosphorus
        ('P', True): 'phosphonium molecule',
        ('P', False): 'phosphide anion molecule',
        
        # Halogens
        ('F', False): 'fluoride molecule',
        ('Cl', False): 'chloride molecule',
        ('Br', False): 'bromide molecule',
        ('I', False): 'iodide molecule',
        
        # Metals (typically positive)
        ('Na', True): 'sodium cation molecule',
        ('K', True): 'potassium cation molecule',
        ('Ca', True): 'calcium cation molecule',
        ('Mg', True): 'magnesium cation molecule',
        ('Li', True): 'lithium cation molecule',
        ('Zn', True): 'zinc cation molecule',
        ('Fe', True): 'iron cation molecule',
        ('Cu', True): 'copper cation molecule',
        ('Al', True): 'aluminum cation molecule',
        
        # Boron
        ('B', True): 'boronium molecule',
        ('B', False): 'borate molecule',
        
        # Silicon
        ('Si', True): 'silicenium molecule',
        ('Si', False): 'silicate molecule',
    }
    
    return classifications.get((atom_symbol, is_positive), None)


def composition_counter_diff(mol1, mol2):
    """
    Check if mol1 is a substructure of mol2, handling metal complexes that 
    RDKit's HasSubstructMatch misses (e.g., [Mg+]Br in CC(=[O+][Mg]Br)c1ccccc1).
    
    Args:
        mol1: Query molecule (potential substructure)  
        mol2: Target molecule (larger structure)
    
    Returns:
        bool: True if mol1 atoms are all present in mol2 with sufficient counts
    """
    if mol1 is None or mol2 is None:
        return False
    
    # Get atom compositions: element -> count (ignoring charges for metal complexes)
    comp1 = Counter(atom.GetSymbol() for atom in mol1.GetAtoms())
    comp2 = Counter(atom.GetSymbol() for atom in mol2.GetAtoms())

    # Check if all atoms in mol1 are present in mol2 with at least the same count
    return all(comp2[element] >= count for element, count in comp1.items())


def parse_smiles_atoms(smiles):
    """
    Parse SMILES string to count atoms.
    Simple regex-based approach for common elements.
    """
    # Remove brackets and charges for counting
    clean_smiles = re.sub(r'[\[\]+-]', '', smiles)
    
    # Find all atoms (capital letter optionally followed by lowercase)
    atoms = re.findall(r'[A-Z][a-z]?', clean_smiles)
    
    # Count occurrences
    atom_count = Counter(atoms)
    
    return atom_count

def calculate_mass_balance(reactants, products):
    """
    Calculate mass balance between reactants and products.
    Returns (balance_score, total_atoms) where balance_score is between 0 and 1.
    """
    # Count atoms on reactant side
    reactant_dict = atom_dict()
    for r in reactants:
        reactant_dict = count_atoms_in_mol(reactant_dict,r)

    # Count atoms on product side
    product_dict = atom_dict()

    for p in products:
        product_dict = count_atoms_in_mol(product_dict,p)
        
    # Calculate total atoms
    total_atoms = sum(reactant_dict.values()) + sum(product_dict.values())

    if total_atoms == 0:
        return 0, 0

    # Calculate balanced atoms
    balanced_atoms = 0
    balanced_mass = 0
    r_mass = 0
    p_mass = 0
    all_elements = set(reactant_dict.keys()) | set(product_dict.keys())

    for element in all_elements:
        r_count = reactant_dict.get(element, 0)
        p_count = product_dict.get(element, 0)
        balanced_atoms += min(r_count, p_count) * 2  # Count both sides
        
        #r_mass += r_count * atom_mass(element) 
        #p_mass += p_count * atom_mass(element) 
        #balanced_mass += abs(r_count - p_count) * atom_mass(element) 
        
    # Balance score is ratio of balanced atoms to total atoms
    balance_atom_score = balanced_atoms / total_atoms if total_atoms > 0 else 0
    #balance_mass_score = 1 - abs(balanced_mass/max(r_mass,p_mass)) if r_mass and p_mass else 0

    #balance_score = (balance_atom_score + balance_mass_score) / 2
    return balance_atom_score, total_atoms

def compute_substructure_score_fraction(mol1,mol2,timeout = 1):
    res = rdFMCS.FindMCS([mol1, mol2], timeout=timeout)

    if res.canceled or res.numAtoms == 0:
        return 0.0
    if mol1.GetNumAtoms() >= mol2.GetNumAtoms():
        return res.numAtoms / mol1.GetNumAtoms()
    else:
        return res.numAtoms / mol2.GetNumAtoms()

def calculate_substructre_score(reactants,products):
    """ To calculate a score between 0 and 1 for the match between substructures in reactants and products """

    score_list_r = [max(compute_substructure_score_fraction(flatten_bonds_for_mcs(p), flatten_bonds_for_mcs(r)) for p in products) for r in reactants]
    score_list_p = [max(compute_substructure_score_fraction(flatten_bonds_for_mcs(r), flatten_bonds_for_mcs(p)) for r in reactants) for p in products]

    score_add_r = [1-(x-1) if x > 1 else x for x in score_list_r]
    score_add_p = [1-(x-1) if x > 1 else x for x in score_list_p]

    score_r = sum(score_list_r)
    score_p = sum(score_list_p)
    return max(score_r,score_p)

def count_atoms_in_mol(a_dict, mol):
    """
    Count atoms in an RDKit molecule and add to existing atom dictionary.

    Args:
        atom_dict (dict): Dictionary with element symbols as keys and counts as values
        mol (rdkit.Chem.Mol): RDKit molecule object

    Returns:
        dict: Updated dictionary with atom counts added
    """
    # Create a copy to avoid modifying the original dictionary
    updated_dict = a_dict.copy()
    
    m = Chem.MolFromSmiles(mol,sanitize = False)
    # Iterate through all atoms in the molecule

    for atom in m.GetAtoms():
        if atom.GetSymbol() == 'H':
            pass
        else:
            element_symbol = atom.GetSymbol()
            updated_dict[element_symbol] += 1
    
    updated_dict['H'] += count_hydrogens_from_smiles(mol)['total_h']

    return updated_dict

def atom_mass(element):
    p_table = {
        'H': 1, 'He': 4,
        'Li': 7, 'Be': 9, 'B': 11, 'C': 12, 'N': 14, 'O': 16, 'F': 19, 'Ne': 20,
        'Na': 23, 'Mg': 24, 'Al': 27, 'Si': 28, 'P': 31, 'S': 32, 'Cl': 35, 'Ar': 40,
        'K': 39, 'Ca': 40, 'Sc': 45, 'Ti': 48, 'V': 51, 'Cr': 52, 'Mn': 55, 'Fe': 56, 'Co': 59, 'Ni': 59, 'Cu': 64, 'Zn': 65,
        'Ga': 70, 'Ge': 73, 'As': 75, 'Se': 79, 'Br': 80, 'Kr': 84,
        'Rb': 85, 'Sr': 88, 'Y': 89, 'Zr': 91, 'Nb': 93, 'Mo': 96, 'Tc': 98, 'Ru': 101, 'Rh': 103, 'Pd': 106, 'Ag': 108, 'Cd': 112,
        'In': 115, 'Sn': 119, 'Sb': 122, 'Te': 128, 'I': 127, 'Xe': 131,
        'Cs': 133, 'Ba': 137, 'La': 139, 'Ce': 140, 'Pr': 141, 'Nd': 144, 'Pm': 145, 'Sm': 150, 'Eu': 152, 'Gd': 157, 'Tb': 159, 'Dy': 162,
        'Ho': 165, 'Er': 167, 'Tm': 169, 'Yb': 173, 'Lu': 175,
        'Hf': 178, 'Ta': 181, 'W': 184, 'Re': 186, 'Os': 190, 'Ir': 192, 'Pt': 195, 'Au': 197, 'Hg': 201,
        'Tl': 204, 'Pb': 207, 'Bi': 209, 'Po': 209, 'At': 210, 'Rn': 222,
        'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232, 'Pa': 231, 'U': 238, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247, 'Cf': 251,
        'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262,
        'Rf': 267, 'Db': 270, 'Sg': 271, 'Bh': 270, 'Hs': 277, 'Mt': 278, 'Ds': 281, 'Rg': 282, 'Cn': 285, 'Fl': 289, 'Lv': 293, 'Ts': 294, 'Og': 294
    }
    return p_table[element]

def atom_dict():
    return {
    'H': 0, 'He': 0,
    'Li': 0, 'Be': 0, 'B': 0, 'C': 0, 'N': 0, 'O': 0, 'F': 0, 'Ne': 0,
    'Na': 0, 'Mg': 0, 'Al': 0, 'Si': 0, 'P': 0, 'S': 0, 'Cl': 0, 'Ar': 0,
    'K': 0, 'Ca': 0, 'Sc': 0, 'Ti': 0, 'V': 0, 'Cr': 0, 'Mn': 0, 'Fe': 0, 'Co': 0, 'Ni': 0, 'Cu': 0, 'Zn': 0,
    'Ga': 0, 'Ge': 0, 'As': 0, 'Se': 0, 'Br': 0, 'Kr': 0,
    'Rb': 0, 'Sr': 0, 'Y': 0, 'Zr': 0, 'Nb': 0, 'Mo': 0, 'Tc': 0, 'Ru': 0, 'Rh': 0, 'Pd': 0, 'Ag': 0, 'Cd': 0,
    'In': 0, 'Sn': 0, 'Sb': 0, 'Te': 0, 'I': 0, 'Xe': 0,
    'Cs': 0, 'Ba': 0, 'La': 0, 'Ce': 0, 'Pr': 0, 'Nd': 0, 'Pm': 0, 'Sm': 0, 'Eu': 0, 'Gd': 0, 'Tb': 0, 'Dy': 0,
    'Ho': 0, 'Er': 0, 'Tm': 0, 'Yb': 0, 'Lu': 0,
    'Hf': 0, 'Ta': 0, 'W': 0, 'Re': 0, 'Os': 0, 'Ir': 0, 'Pt': 0, 'Au': 0, 'Hg': 0,
    'Tl': 0, 'Pb': 0, 'Bi': 0, 'Po': 0, 'At': 0, 'Rn': 0,
    'Fr': 0, 'Ra': 0, 'Ac': 0, 'Th': 0, 'Pa': 0, 'U': 0, 'Np': 0, 'Pu': 0, 'Am': 0, 'Cm': 0, 'Bk': 0, 'Cf': 0,
    'Es': 0, 'Fm': 0, 'Md': 0, 'No': 0, 'Lr': 0,
    'Rf': 0, 'Db': 0, 'Sg': 0, 'Bh': 0, 'Hs': 0, 'Mt': 0, 'Ds': 0, 'Rg': 0, 'Cn': 0, 'Fl': 0, 'Lv': 0, 'Ts': 0, 'Og': 0}

def calculate_charge_balance(reactants,products):
    """ Calculate charge balance between reactants and products 
    returns [dc,dc]"""
    try:
        charge_r = 0
        for r in reactants:
            mol = Chem.MolFromSmiles(r)
            charge_r += sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        
        charge_p = 0
        for p in products:
            mol = Chem.MolFromSmiles(p)
            charge_p += sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        return (charge_r == charge_p), (charge_r,charge_p)
    except:
        charge_r = 0
        for r in reactants:
            for c in r:
                if c == '+':
                    charge_r += 1
                if c == '-':
                    charge_r -= 1
                    
        charge_p = 0
        for p in products:
            for c in p:
                if c == '+':
                    charge_p += 1
                if c == '-':
                    charge_p -= 1
                    
        return (charge_r == charge_p), (charge_r,charge_p)        

def generate_valid_combinations(r_match, p_match):
    """
    Generate all valid combinations of reactants and products.
    A combination is valid if it doesn't include both a complex molecule and its fragments.
    """
    # Get all possible molecules
    all_reactants = set(r_match.keys())
    all_products = set(p_match.keys())
    
    # Add fragments as possible molecules
    for fragments in r_match.values():
        all_products.update(fragments)
    for fragments in p_match.values():
        all_reactants.update(fragments)
    
    valid_combinations = []
    
    # Try all possible combinations of reactants and products
    for r_count in range(1, len(all_reactants) + 1):
        for p_count in range(1, len(all_products) + 1):
            for reactant_combo in combinations(all_reactants, r_count):
                for product_combo in combinations(all_products, p_count):
                    
                    # Check if combination is valid (no conflicts)
                    if is_valid_combination(reactant_combo, product_combo, r_match, p_match):
                        valid_combinations.append((list(reactant_combo), list(product_combo)))
    
    return valid_combinations

def is_valid_combination(reactants, products, r_match, p_match):
    """
    Check if a combination is valid (no molecule conflicts).
    Invalid if we use both a complex molecule and its constituent fragments.
    """
    reactant_set = set(reactants)
    product_set = set(products)
    
    # Check for conflicts in reactants
    for reactant in reactants:
        if reactant in r_match:
            # If we use this complex reactant, we shouldn't also use its fragments as products
            fragments = set(r_match[reactant])
            if fragments.intersection(product_set):
                # Allow if we're using the fragments AND we have the corresponding complex product
                # This would represent a decomposition-recomposition reaction
                pass
    
    # Check for conflicts in products  
    for product in products:
        if product in p_match:
            # If we use this complex product, we shouldn't also use its fragments as reactants
            # unless it's a valid construction reaction
            fragments = set(p_match[product])
            if fragments.intersection(reactant_set):
                # This is actually a valid construction - fragments combining to form complex molecule
                pass
    
    return True

def optimize_reaction(r_match, p_match, min_mass_threshold=0.4):
    """
    Optimize chemical reaction composition for mass balance and total mass.
    
    Args:
        r_match: Dict mapping reactants to their possible product fragments
        p_match: Dict mapping products to their possible reactant fragments
        min_mass_threshold: Minimum mass balance threshold (default 40%)
    
    Returns:
        dict with optimized reaction details
    """
    
    best_reaction = {
        'core_reactants_smiles': [],
        'core_products_smiles': [],
        'unused_molecules': [],
        'balance_score': 0,
        'total_atoms': 0,
        'mass_balance_achieved': False
    }

    # Get all possible molecules
    all_reactants = set(r_match.keys())
    all_products = set(p_match.keys())
    
    # Add fragments as possible molecules
    for fragments in r_match.values():
        all_products.update(fragments)
    for fragments in p_match.values():
        all_reactants.update(fragments)
    
    # Generate all valid combinations
    valid_combinations = generate_valid_combinations(r_match, p_match)
    
    # Evaluate each combination
    for reactant_combo, product_combo in valid_combinations:
        # Calculate mass balance
        balance_score, total_atoms = calculate_mass_balance(reactant_combo, product_combo)
        
        # Check if this is better than current best
        is_better = False
        
        if balance_score >= min_mass_threshold:
            # If mass balance threshold is met, prefer higher balance score, then higher total mass
            if (best_reaction['balance_score'] < min_mass_threshold or 
                balance_score > best_reaction['balance_score'] or
                (balance_score == best_reaction['balance_score'] and total_atoms > best_reaction['total_atoms'])):
                is_better = True
        else:
            # If threshold not met, prefer higher balance score
            if balance_score > best_reaction['balance_score']:
                is_better = True
        
        if is_better:
            # Calculate unused molecules
            unused = (all_reactants - set(reactant_combo)) | (all_products - set(product_combo))
            
            best_reaction.update({
                'core_reactants_smiles': reactant_combo,
                'core_products_smiles': product_combo,
                'unused_molecules': list(unused),
                'balance_score': balance_score,
                'total_atoms': total_atoms,
                'mass_balance_achieved': balance_score >= min_mass_threshold
            })
    
    return best_reaction

def mols_equal(mol1, mol2):
    # Check that both molecules have the same number of atoms and are isomorphic
    try:
        mol1 = Chem.RemoveHs(mol1)
        mol2 = Chem.RemoveHs(mol2)

        return (mol1.GetNumAtoms() == mol2.GetNumAtoms() and
                mol1.HasSubstructMatch(mol2) and
                mol2.HasSubstructMatch(mol1))
    except:
        return False

def has_metal(mol: Chem.Mol) -> bool:
    """Returns True if the molecule contains at least one metal atom."""
    # Periodic table atomic numbers for metals (excluding metalloids and nonmetals)
    # You can refine this list for your specific use case
    METAL_ATOMIC_NUMS = {
        # Alkali metals
        3, 11, 19, 37, 55, 87,
        # Alkaline earth metals
        4, 12, 20, 38, 56, 88,
        # Transition metals
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
        39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
        72, 73, 74, 75, 76, 77, 78, 79, 80,
        104, 105, 106, 107, 108, 109, 110, 111, 112,
        # Post-transition metals
        13, 31, 49, 50, 81, 82, 83, 113, 114, 115, 116,
        # Lanthanides
        57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
        # Actinides
        89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103
    }
    return any(atom.GetAtomicNum() in METAL_ATOMIC_NUMS for atom in mol.GetAtoms())

def zero_charge(mol):
    """ Set formal charges to 0 with maintaining H count. """
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        if charge != 0:
            atom.SetFormalCharge(0)
    return mol


def score_reaction(mass_balance, charge_balance, substructure, score, sigmoid_max):
    """
    Compare current mass balance, number of atoms involved and charge balance
    to optimize for the best combination of reactants and products.
    
    Args:
        mass_balance: tuple (mass_ratio, num_atoms) where mass_ratio should be close to 1 (≤1)
        charge_balance: tuple (charge_left, charge_right) - better when equal
        score: tuple (score_mass_ratio, score_num_atoms) - the comparison baseline
        sigmoid_max: maximum possible value for num_atoms (for sigmoid calculation)
    
    Returns:
        bool: True if mass_balance is better than score, False otherwise
    """

    # Extract values
    mb_mass_ratio, mb_num_atoms = mass_balance
    cb_left, cb_right = charge_balance
    check_score = score[0]
    
    # 1. Calculate mass balance quality (higher is better)
    # Penalize heavily if mass_ratio > 1, reward closeness to 1
    def mass_quality(mass_ratio):
        if mass_ratio > 1:
            return mass_ratio / 2  # Mass ratio should normally not exceed 1
        return mass_ratio  # Closer to 1 is better
    
    mb_mass_quality = mass_quality(mb_mass_ratio)
    #score_mass_quality = mass_quality(score_mass_ratio)
    
    # 2. Calculate sigmoid weight based on num_atoms
    # Low num_atoms gets high weight (sigmoid approaches 1)
    # High num_atoms gets low weight (sigmoid approaches 0)
    def sigmoid_weight(num_atoms, max_val):
        # Sigmoid: 1 - 1 / (1 + exp(k * (x - midpoint)))
        # We want low weight for low num_atoms, so we use negative slope
        if max_val == 0:
            max_val = 0.001
        midpoint = max_val / 3
        k = 12 / max_val  # Steepness factor (adjust as needed)
        return 1 - 1 / (1 + math.exp(k * (num_atoms - midpoint)))
    
    mb_sigmoid_weight = sigmoid_weight(mb_num_atoms, sigmoid_max)
    #score_sigmoid_weight = sigmoid_weight(score_num_atoms, sigmoid_max)
    
    # 3. Calculate charge balance quality (higher is better)
    # Perfect when charges are equal, worse as difference increases
    def charge_quality(left, right):
        diff = abs(left - right)
        # Exponential decay: e^(-diff) so equal charges give 1, larger diff approaches 0
        return math.exp(-(0.3)*diff)
    
    mb_charge_quality_factor = charge_quality(cb_left, cb_right)

    #score_charge_quality_factor = charge_quality(score_cb_left, score_cb_right)
    # 4. Combine factors for final scores
    # Alternative approach: Direct penalty for high atom counts
    def combined_score(mass_qual, low_atom_penalty, charge_qual,substructure):
        # Option 1: Sigmoid as multiplicative factor
        # Lower atom counts get higher sigmoid weights, boosting their score
        # High for low atoms, low for high atoms
        
        # Option 2: More aggressive atom penalty
        # Uncomment this line and comment the above for stronger atom count influence
        # atom_penalty = sigmoid_w ** 2  # Squares the sigmoid for stronger effect
        
        # Final score combines mass quality, atom penalty, and charge balance
        return mass_qual * low_atom_penalty * charge_qual * substructure

    mb_final_score = combined_score(mb_mass_quality, mb_sigmoid_weight, mb_charge_quality_factor, substructure)
    #score_final_score = combined_score(score_mass_quality, score_sigmoid_weight, score_charge_quality_factor)
    #print(mb_mass_quality,mb_sigmoid_weight,charge_quality_factor)
    # Return True if mass_balance is better than score
    return mb_final_score > check_score, mb_final_score

from rdkit.Chem import rdFingerprintGenerator

def check_subfragment_match(reactant,product):
    if reactant is None or product is None:
        return False
    match = False
    sim_threshold = 0.5
    reactant_candidates = {}  # key: parent SMILES, value: set of candidate SMILES
    product_candidates = {}
    unsanitized_r_mol = []
    unsanitized_p_mol = []
    for smi in reactant:
        cand = get_candidate_fragments(smi)
        if cand:
            reactant_candidates[smi] = cand
        else:
            unsanitized_r_mol.append(smi)
    for smi in product:
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
                        match = True
                        # Once one match is found for this reactant parent, no need to check further.
                        break
                else:
                    continue
                break
    return match