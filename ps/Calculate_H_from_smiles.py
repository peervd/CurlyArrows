import re
from rdkit import Chem

def count_hydrogens_from_smiles(smiles):
    """
    Calculate hydrogens in a SMILES string with improved parsing.
    
    Args:
        smiles (str): SMILES string
        
    Returns:
        dict: Dictionary with hydrogen count breakdown
    """
    results = {}
    
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return {'total_h': 0, 'explicit_h': 0, 'implicit_h': 0}
        
        explicit_h = 0
        implicit_h = 0
        
        # Count explicit and implicit hydrogens
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()

            if symbol == 'H':
                explicit_h += 1
            else:
                # Calculate implicit hydrogens for non-H atoms
                implicit_h += calculate_implicit_h_for_atom(atom, mol)
        
        results = {
            'explicit_h': explicit_h,
            'implicit_h': implicit_h,
            'total_h': explicit_h + implicit_h
        }
    
    except Exception as e:
        # Fallback to string parsing if RDKit fails
        results = parse_smiles_string_for_hydrogens(smiles)

    return results

def calculate_implicit_h_for_atom(atom, mol):
    """
    Calculate implicit hydrogens for a single atom based on valency rules.
    """
    symbol = atom.GetSymbol()
    
    # Skip hydrogen atoms
    if symbol == 'H':
        return 0
    
    # Get atom properties
    formal_charge = atom.GetFormalCharge()
    degree = atom.GetDegree()  # Number of bonds
    is_aromatic = atom.GetIsAromatic()
    
    # Define valency rules
    valency_rules = {
        'C': {1:[3],0:[4],-1:[3]},
        'N': {1:[4],0:[3],-1:[2]},
        'O': {1:[3],0:[2],-1:[1]},
        'S': {1:[3],0:[2, 4, 6],-1:[1]},
        'P': {1:[4],0:[3, 5],-1:[2]},
        'F': {0:[1],-1:[0]},
        'Cl': {0:[1],-1:[0]},
        'Br': {0:[1],-1:[0]},
        'I': {0:[1],-1:[0]}
    }
    
    if symbol not in valency_rules:
        return 0
    
    # Get possible valencies for this atom
    possible_valencies = valency_rules[symbol]
    
    # Count total bond order (sum of bond orders)
    total_bond_order = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
    
    # For aromatic atoms, adjust calculation
    if is_aromatic:
        # For aromatic carbons, typical valency is 4
        if symbol == 'C':
            expected_valency = 4
        elif symbol == 'N':
            expected_valency = 3
        else:
            expected_valency = possible_valencies[0]

        # Calculate implicit hydrogens
        implicit_h = max(0, expected_valency - total_bond_order - abs(formal_charge))
        return int(implicit_h)
    
    else:
        # For non-aromatic atoms
        # Special handling for charged atoms
        if formal_charge != 0:
            # For positively charged atoms, they need more bonds to satisfy valency
            # For negatively charged atoms, they need fewer bonds
            
            # Find the appropriate valency considering the charge
            for charge in possible_valencies.keys():
                # For positive charge: atom needs more electrons (more bonds/hydrogens)
                # For negative charge: atom has extra electrons (fewer bonds/hydrogens)
                try:
                    if formal_charge == charge:
                        required_bonds = possible_valencies[charge]
                    implicit_h = []
                    for x in required_bonds:
                        if x >= total_bond_order and x - total_bond_order >= 0:
                            implicit_h.append(x - total_bond_order)
                        return(int(max(implicit_h)))
                except:
                    return 0
                    
                
        else:
            # Neutral atom - find the lowest valency that satisfies current bonding
            for valency in possible_valencies[0]:
                if valency >= total_bond_order:
                    implicit_h = max(0, valency - total_bond_order)
                    return int(implicit_h)
    
    return 0

def parse_smiles_string_for_hydrogens(smiles):
    """
    Fallback method: parse SMILES string directly for hydrogen information.
    """
    explicit_h = 0
    
    # Find all bracketed atoms
    bracket_pattern = r'\[([^\]]+)\]'
    matches = re.findall(bracket_pattern, smiles)
    
    for match in matches:
        # Look for H followed by optional number
        h_matches = re.findall(r'H(\d*)', match)
        for h_match in h_matches:
            count = int(h_match) if h_match else 1
            explicit_h += count
    
    # Count standalone H atoms (not in brackets)
    # This is tricky and approximate
    non_bracket_h = 0
    temp_smiles = re.sub(r'\[[^\]]+\]', '', smiles)  # Remove bracketed parts
    non_bracket_h = temp_smiles.count('H')
    
    total_explicit = explicit_h + non_bracket_h
    
    return {
        'explicit_h': total_explicit,
        'implicit_h': 0,  # Can't calculate without proper parsing
        'total_h': total_explicit
    }

def count_atoms_with_improved_h(a_dict, smiles):
    """
    Enhanced atom counting with improved hydrogen calculation.
    """
    updated_dict = a_dict.copy()
    
    # Get hydrogen count
    h_results = count_hydrogens_from_smiles(smiles)
    
    # Count all atoms
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                if symbol not in updated_dict:
                    updated_dict[symbol] = 0
                updated_dict[symbol] += 1
    except:
        pass
    
    # Set hydrogen count
    if 'H' not in updated_dict:
        updated_dict['H'] = 0
    updated_dict['H'] = h_results['total_h']
    
    return updated_dict
