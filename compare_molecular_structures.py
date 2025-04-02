from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

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
            return ["does not match the model answer_0",carbon_count]
            
    except Exception as e:
        return f"Error comparing molecules: {str(e)}"