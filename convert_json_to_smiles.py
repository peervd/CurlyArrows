import json
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, rdMolAlign, rdmolops
import matplotlib.pyplot as plt

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def getSubmolRadAtom(mol, index, radius):
    """
    Extract a submol around a specific atom with a given radius
    
    Parameters:
    mol (Chem.Mol): RDKit molecule
    index (int): Atom index to center the submol extraction
    radius (int): Radius of atom environment to extract
    
    Returns:
    str: SMILES representation of the extracted submol
    """
    # Find atom environment
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, index)
    
    # Prepare atom map for submol extraction
    amap = {}
    
    # Extract submol
    submol = Chem.PathToSubmol(mol, env, atomMap=amap)
    
    # If no submol found, try with smaller radius
    if not bool(amap):
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 1, index)
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
    
    # Convert submol to SMILES
    subsmi = Chem.MolToSmiles(submol)
    return subsmi

def getSubmolRadBond(mol, index, radius):
    """
    Extract a submol around a specific bond with a given radius
    
    Parameters:
    mol (Chem.Mol): RDKit molecule
    index (tuple): Tuple of atom indices defining the bond
    radius (int): Radius of bond environment to extract
    
    Returns:
    str: SMILES representation of the extracted submol
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
    
    # Convert submol to SMILES
    submol_smiles = Chem.MolToSmiles(submol)
    
    return submol_smiles

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

def parse_json_to_smiles(json_string, struc = False):
    data = json.loads(json_string)
    molecules = data.get("m", [])
    synthetic_arrows = sorted([(arrow["x1"], arrow["x2"], (arrow["y1"] + arrow["y2"]) / 2) for arrow in data.get("s", []) if arrow["t"] == "Line" and arrow.get("a") == "synthetic"],key=lambda x: x[0])
    resonance_arrows = sorted([(arrow["y1"], arrow["y2"], (arrow["x1"] + arrow["x2"]) / 2) for arrow in data.get("s", []) if arrow["t"] == "Line" and arrow.get("a") == "resonance"],key=lambda x: x[0])
    curved_arrows = [arrow for arrow in data.get("s", []) if arrow["t"] == "Pusher"]
    smiles_dict = {"reactants": {"molecules": {}, "arrows": []}, "products": {"molecules": {}, "arrows": []}}
    mol_struc = {"reactants": {}, "products": {}}

    for i in range(len(synthetic_arrows) - 1):
        if any(synthetic_arrows[i][1] < x[2] < synthetic_arrows[i+1][0] for x in resonance_arrows):
            smiles_dict[f"intermediates_{i+1}_a"] = {"molecules": {}, "arrows": []}
            mol_struc[f"intermediates_{i+1}_a"] = {}
            smiles_dict[f"intermediates_{i+1}_b"] = {"molecules": {}, "arrows": []}
            mol_struc[f"intermediates_{i+1}_b"] = {}             
        else:    
            smiles_dict[f"intermediates_{i+1}"] = {"molecules": {}, "arrows": []}
            mol_struc[f"intermediates_{i+1}"] = {} 

    molecule_categories = {}
    atom_map = {}
    atom_ids = {}
    bond_ids = {}    
    m = {}    
    for mol_number, mol_data in enumerate(molecules):
        atom_map[mol_number] = {}
        atom_ids[mol_number] = {}
        bond_ids[mol_number] = {}    
        
        mol = Chem.RWMol()
        
        # Determine molecule position
        x_positions = [atom["x"] for atom in mol_data.get("a", [])]
        y_positions = [atom["y"] for atom in mol_data.get("a", [])]
        avg_x = sum(x_positions) / len(x_positions) if x_positions else 0
        avg_y = sum(y_positions) / len(y_positions) if y_positions else 0
        
        # Add atoms with formal charges and explicit hydrogens
        for idx, atom in enumerate(mol_data.get("a", [])):
            element = atom.get("l", "C")  # Default to carbon if no label
            rd_atom = Chem.Atom(element)
            if "c" in atom:
                rd_atom.SetFormalCharge(atom["c"])  # Set formal charge
            atom_idx = mol.AddAtom(rd_atom)
            atom_map[mol_number][idx] = atom_idx
            if "i" in atom:
                atom_ids[mol_number][atom["i"]] = atom_idx

        # Compute implicit valence before adding bonds
        mol.UpdatePropertyCache(strict=False)
        m_seq = Chem.MolToSmiles(mol)

        # Add bonds
        for bond in mol_data.get("b", []):
            start = atom_map[mol_number][bond["b"]]
            end = atom_map[mol_number][bond["e"]]
            order = bond.get("o", 1)  # Default single bond
            
            if order == 1:
                bond_type = rdchem.BondType.SINGLE
            elif order == 2:
                bond_type = rdchem.BondType.DOUBLE
            elif order == 3:
                bond_type = rdchem.BondType.TRIPLE
            else:
                continue
            
            mol.AddBond(start, end, bond_type)

            if "i" in bond:
                bond_ids[mol_number][bond["i"]] = (start, end)
        
        
        # Generate SMILES with explicit hydrogen and charges
        sanitization = True
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            sanitization = False
            #print("Sanitization failed:", e)
        m[mol_number] = mol
        smiles = Chem.MolToSmiles(mol)

        # Categorize molecules based on position relative to arrows
        if not synthetic_arrows or avg_x < synthetic_arrows[0][0]:
            category = "reactants"
        elif avg_x > synthetic_arrows[-1][1]:
            category = "products"
        else:
            for i, (x1, x2, y) in enumerate(synthetic_arrows[:-1]):
                if any(x2 < x[2] < synthetic_arrows[i+1][0] for x in resonance_arrows):
                    if x2 < avg_x < synthetic_arrows[i+1][0] and any(avg_y < y[0] for y in resonance_arrows):
                        category = f"intermediates_{i+1}_a"
                        break
                    elif x2 < avg_x < synthetic_arrows[i+1][0] and any(avg_y > y[1] for y in resonance_arrows):
                        category = f"intermediates_{i+1}_b"
                        break
                elif x2 < avg_x < synthetic_arrows[i+1][0]:
                    category = f"intermediates_{i+1}"
                    break

        smiles_dict[category]["molecules"][mol_number] = smiles
        molecule_categories[mol_data["a"][0]["i"]] = category
        mol_struc[category][mol_number] = sanitization

    # Process curved arrows
    for arrow in curved_arrows:
        start = arrow['o1']
        end = arrow['o2']
        start_a = None
        end_a = None
        for x in atom_ids:
            if start in atom_ids[x]:
                mol_start = x
                start_a = atom_ids[x].get(start)
                len_start = len(atom_ids[x])
                s = 'atom'
            if end in atom_ids[x]:
                mol_end = x
                end_a = atom_ids[x].get(end)
                len_end = len(atom_ids[x])
                e = 'atom'
        for x in bond_ids:
            if start in bond_ids[x]:
                mol_start = x
                start_a = bond_ids[x].get(start)
                len_start = len(atom_ids[x])
                s = 'bond'
            if end in bond_ids[x]:
                mol_end = x
                end_a = bond_ids[x].get(end)
                len_end = len(atom_ids[x])
                e = 'bond'
           
        SubMol = {}
        if s == 'atom':
            start_arrow = m[mol_start].GetAtomWithIdx(start_a).GetSymbol()
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadAtom(m[mol_start], start_a, 2)
        elif s == 'bond':
            
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
                start_arrow = Chem.MolToSmiles(m[mol_start])
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadBond(m[mol_start], start_a, 2)
                start_arrow = getRadBond(m[mol_start], start_a)


        if e == 'atom':
            end_arrow = m[mol_end].GetAtomWithIdx(end_a).GetSymbol()
            if len_end <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
            elif len_end > 2:
                SubMol['end_sub'] = getSubmolRadAtom(m[mol_end], end_a, 2)            
        elif e == 'bond':
            if len_start <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
                end_arrow = Chem.MolToSmiles(m[mol_end])
            elif len_start > 2:
                SubMol['end_sub'] = getSubmolRadBond(m[mol_end], end_a, 2)
                end_arrow = getRadBond(m[mol_end], end_a)

        for key in smiles_dict:
            if mol_start in smiles_dict[key]['molecules'].keys():
                set_key = key
        if start is not None and end is not None:
            smiles_dict[set_key]["arrows"].append({
                "start_molecule": mol_start,
                "start_arrow": start_arrow,
                "substructure_start":  SubMol['start_sub'], 
                "end_molecule": mol_end,
                "end_arrow": end_arrow,
                "substructure_end": SubMol['end_sub'],
                "electrons": arrow["e"]
            })
    
    if struc == True:
        return smiles_dict,mol_struc
    elif struc == False:
        return smiles_dict

