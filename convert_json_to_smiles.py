import json
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem, rdMolAlign, rdmolops
import matplotlib.pyplot as plt

def getSubmolRadAtom(mol, index, radius):
    env=Chem.FindAtomEnvironmentOfRadiusN(mol, radius, index)
    amap={}
    submol=Chem.PathToSubmol(mol, env, atomMap=amap)
    if bool(amap) == False:
        env=Chem.FindAtomEnvironmentOfRadiusN(mol, 1, index)
        submol=Chem.PathToSubmol(mol, env, atomMap=amap)
    subsmi=Chem.MolToSmiles(submol)
    return subsmi

def getSubmolRadBond(mol, index, radius):
    ''''retrieve indices of atoms from substructure in molecule'''
    m = []
    amap={}
    for a,ind in enumerate(index):
        env=Chem.FindAtomEnvironmentOfRadiusN(mol, 1, ind)
        submol=Chem.PathToSubmol(mol, env, atomMap=amap)
        m.append(submol)

    matches = mol.GetSubstructMatches(m[0])
    m_list = []
    for x in matches[0]:
        m_list.append(x)
    if len(matches) > 1:
        for x in matches[1]:
            if x not in m_list:
                m_list.append(x)
    m_list.sort()
    mol = Chem.Mol(mol)
    submol = Chem.MolFragmentToSmiles(mol, m_list)
    return submol

def getRadBond(mol, index):
    submol = Chem.MolFragmentToSmiles(mol, index)
    return submol

def parse_json_to_smiles(json_string):
    data = json.loads(json_string)
    molecules = data.get("m", [])
    arrows = sorted([(arrow["x1"], arrow["x2"]) for arrow in data.get("s", []) if arrow["t"] == "Line"], key=lambda x: x[0])
    curved_arrows = [arrow for arrow in data.get("s", []) if arrow["t"] == "Pusher"]
    
    smiles_dict = {"reactants": {"molecules": {}, "arrows": []}, "products": {"molecules": {}, "arrows": []}}
    
    for i in range(len(arrows) - 1):
        smiles_dict[f"intermediates_{i+1}"] = {"molecules": {}, "arrows": []}
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
        avg_x = sum(x_positions) / len(x_positions) if x_positions else 0
        
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
        Chem.SanitizeMol(mol)
        m[mol_number] = mol
        smiles = Chem.MolToSmiles(mol)

        # Categorize molecules based on position relative to arrows
        if not arrows or avg_x < arrows[0][0]:
            category = "reactants"
        elif avg_x > arrows[-1][1]:
            category = "products"
        else:
            for i, (x1, x2) in enumerate(arrows[:-1]):
                if x1 < avg_x < arrows[i+1][0]:
                    category = f"intermediates_{i+1}"
                    break
        
        smiles_dict[category]["molecules"][mol_number] = smiles
        molecule_categories[mol_data["a"][0]["i"]] = category

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
            start_arrow = m[mol_start].GetAtomWithIdx(start_a).GetSymbol(),start_a
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadAtom(m[mol_start], start_a, 2)
        elif s == 'bond':
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
                start_arrow = Chem.MolToSmiles(m[mol_start]),start_a
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadBond(m[mol_start], start_a, 2)
                start_arrow = getRadBond(m[mol_start], start_a),start_a
            
        if e == 'atom':
            end_arrow = m[mol_end].GetAtomWithIdx(end_a).GetSymbol(),end_a
            if len_end <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
            elif len_end > 2:
                SubMol['end_sub'] = getSubmolRadAtom(m[mol_end], end_a, 2)            
        elif e == 'bond':
            if len_start <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
                end_arrow = Chem.MolToSmiles(m[mol_end]),end_a
            elif len_start > 2:
                SubMol['end_sub'] = getSubmolRadBond(m[mol_end], end_a, 2)
                end_arrow = getRadBond(m[mol_end], end_a),end_a

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
            
    
    return smiles_dict

def parse_json_to_reaction_rule(name,json_string):
    data = json.loads(json_string)
    molecules = data.get("m", [])
    arrows = sorted([(arrow["x1"], arrow["x2"]) for arrow in data.get("s", []) if arrow["t"] == "Line"], key=lambda x: x[0])
    curved_arrows = [arrow for arrow in data.get("s", []) if arrow["t"] == "Pusher"]
    string = [name]
    smiles_dict = {"reactants": {"molecules": {}, "arrows": []}, "products": {"molecules": {}, "arrows": []}}
    
    for i in range(len(arrows) - 1):
        smiles_dict[f"intermediates_{i+1}"] = {"molecules": {}, "arrows": []}
    molecule_categories = {}
    atom_map = {}
    atom_ids = {}
    bond_ids = {}    
    m = {}   
    c_r = 0 
    c_p = 0
    for mol_number, mol_data in enumerate(molecules):
        atom_map[mol_number] = {}
        atom_ids[mol_number] = {}
        bond_ids[mol_number] = {}    
        
        mol = Chem.RWMol()
        
        # Determine molecule position
        x_positions = [atom["x"] for atom in mol_data.get("a", [])]
        avg_x = sum(x_positions) / len(x_positions) if x_positions else 0
        
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
        Chem.SanitizeMol(mol)
        m[mol_number] = mol
        smiles = Chem.MolToSmiles(mol)
        
        # Categorize molecules based on position relative to arrows
        if not arrows or avg_x < arrows[0][0]:
            category = "reactants"
            c =+ 1
        elif avg_x > arrows[-1][1]:
            category = "products"
        else:
            for i, (x1, x2) in enumerate(arrows[:-1]):
                if x1 < avg_x < arrows[i+1][0]:
                    category = f"intermediates_{i+1}"
                    break
 
        if category == "reactants" and c_r ==0:
            r0 = smiles
            c_r =+ 1
        elif category == "reactants" and c_r > 0:
            r0 = r0+'.'+smiles
        if category == "products" and c_p ==0:
            r1 = smiles
            c_p =+ 1
        elif category == "products" and c_p > 0:
            r1 = r1+'.'+smiles
                
        smiles_dict[category]["molecules"][mol_number] = smiles
        molecule_categories[mol_data["a"][0]["i"]] = category
    string.append(r0)
    string.append(r1)
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
            elif end in atom_ids[x]:
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
            elif end in bond_ids[x]:
                mol_end = x
                end_a = bond_ids[x].get(end)
                len_end = len(atom_ids[x])
                e = 'bond'
        SubMol = {}
        if s == 'atom':
            start_arrow = m[mol_start].GetAtomWithIdx(start_a).GetSymbol(),start_a
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadAtom(m[mol_start], start_a, 2)
        elif s == 'bond':
            if len_start <= 2:
                SubMol['start_sub'] = Chem.MolToSmiles(m[mol_start])
                start_arrow = Chem.MolToSmiles(m[mol_start]),start_a
            elif len_start > 2:
                SubMol['start_sub'] = getSubmolRadBond(m[mol_start], start_a, 2)
                start_arrow = getRadBond(m[mol_start], start_a),start_a
                
        if e == 'atom':
            end_arrow = m[mol_end].GetAtomWithIdx(end_a).GetSymbol(),end_a
            if len_end <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
            elif len_end > 2:
                SubMol['end_sub'] = getSubmolRadAtom(m[mol_end], end_a, 2)            
        elif e == 'bond':
            if len_start <= 2:
                SubMol['end_sub'] = Chem.MolToSmiles(m[mol_end])
                end_arrow = Chem.MolToSmiles(m[mol_end]),end_a
            elif len_start > 2:
                SubMol['end_sub'] = getSubmolRadBond(m[mol_end], end_a, 2)
                end_arrow = getRadBond(m[mol_end], end_a),end_a

        for key in smiles_dict:
            if mol_start in smiles_dict[key]['molecules'].keys():
                set_key = key
        if start is not None and end is not None:
            string.append([mol_start,start_arrow])
            string.append([mol_end,end_arrow])

    return string
