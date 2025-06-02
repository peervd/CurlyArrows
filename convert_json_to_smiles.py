import json
from rdkit import Chem
from rdkit.Chem import rdchem
from molecular_structures import getSubmolRadAtom, getSubmolRadBond, getRadBond
from helper_functions import validate_arrows,sort_list_tuple,pathway_difference

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def molecule_key(molecules, synthetic, resonance, main, index, res_index, category = 'model'):
    if main[index] == True:
        s0 = synthetic[index][1]
        s1 = synthetic[index+2][0]
        y0 = synthetic[index][2]
    elif main[index] == 'pass':
        return 'pass', res_index
    else:
        s0 = synthetic[index][1]
        s1 = synthetic[index+1][0]
        y0 = synthetic[index][2]
    
    # Check if there's a resonance arrow between this synthetic arrow and the next one
    is_resonance = any(s0 < x[2] < s1 for x in resonance)
    if is_resonance and category == 'model':
        return f"intermediates_{index+1}_ra",f"intermediates_{index+1}_rb",(s0,s1,y0), res_index  -1
    if is_resonance and category == 'student':
        return f"intermediates_{index+1}_ra",f"intermediates_{index+1}_rb",(s0,s1,y0), res_index
    
    else:
        # For non-resonance intermediates, check if there's a lower pathway
        has_lower_pathway = False
        # Check if any molecule's average y position is significantly below the pathway
        for mol_data in molecules:
            x_positions = [atom["x"] for atom in mol_data.get("a", [])]
            y_positions = [atom["y"] for atom in mol_data.get("a", [])]
            if not x_positions or not y_positions:
                continue
            avg_x = sum(x_positions) / len(x_positions)
            avg_y = sum(y_positions) / len(y_positions)
            y_arrows = []
            for y_val in synthetic:
                y_arrows.append(y_val[2])
            multiple = pathway_difference(y_arrows)
            # Check if this molecule is between current and next synthetic arrow horizontally
            # and if it's below the average y position of synthetic arrows
            if s0 < avg_x < s1:
                if avg_y > synthetic[index][2] + 50 and multiple == True:  # Threshold for considering "below main chain"
                    has_lower_pathway = True
                    break

        if has_lower_pathway:
            return f"intermediates_{index+1}_ma",f"intermediates_{index+1}_mb", (s0,s1,y0), res_index   
        else:
            return f"intermediates_{res_index+1}",(s0,s1,y0), res_index  
    
def parse_json_to_smiles(json_string, mech_category='model'):
    data = json.loads(json_string)
    molecules = data.get("m", [])
    
    arrows = sorted([(arrow["a"],arrow["x1"], arrow["x2"], arrow["y1"], arrow["y2"]) for arrow in data.get("s", []) if arrow["t"] == "Line" and arrow.get("a") == "synthetic" or arrow.get("a") == "equilibrium" or arrow.get("a") == "resonance" ], key=lambda x: x[0])

    if mech_category == 'student':
        ''' check if the mechanism has the correct structure before further analysis '''
        valid_mech = validate_arrows(arrows)
        if any(x == False for x in valid_mech.values()):
            return 'invalid_mech', valid_mech
    
    reaction_arrows = []
    for a in arrows:
        if a[0] == 'resonance':
            continue
        else:
            reaction_arrows.append((a[0],a[1]))
    
    reaction_arrows = sort_list_tuple(reaction_arrows)

    synthetic_arrows = sorted([(arrow["x1"], arrow["x2"], (arrow["y1"] + arrow["y2"]) / 2) for arrow in data.get("s", []) if arrow["t"] == "Line" and arrow.get("a") == "synthetic" or arrow.get("a") == "equilibrium"], key=lambda x: x[0])
    resonance_arrows = sorted([(arrow["y1"], arrow["y2"], (arrow["x1"] + arrow["x2"]) / 2) for arrow in data.get("s", []) if arrow["t"] == "Line" and arrow.get("a") == "resonance"], key=lambda x: x[0])
    curved_arrows = [arrow for arrow in data.get("s", []) if arrow["t"] == "Pusher"]

    # Initialize dictionary with new naming convention

    smiles_dict = {"reactants": {"molecules": {}, "arrows": []}, "products": {"molecules": {}, "arrows": []}}
    mol_struc = {"reactants": {}, "products": {}}
    
    # check if more than main path is present
    row_syn = []
    double = False
    for i,x in enumerate(synthetic_arrows):
        if i + 1 == len(synthetic_arrows):
            break
        if double == True:
            row_syn.append('pass')
            double = False
            continue
        if synthetic_arrows[i+1][0] < x[0] < synthetic_arrows[i+1][1] or synthetic_arrows[i+1][0] < x[1] < synthetic_arrows[i+1][1]:
            if x[2] * 1.5 < synthetic_arrows[i+1][2]:
                double = True
        elif synthetic_arrows[i+1][0] < x[0] < synthetic_arrows[i+1][1] and synthetic_arrows[i+1][0] < x[1] < synthetic_arrows[i+1][1]:
            if x[2] * 1.5 < synthetic_arrows[i+1][2]:
                double = True
        elif x[0] < synthetic_arrows[i+1][0] and x[1] > synthetic_arrows[i+1][1]:
            if x[2] * 1.5 < synthetic_arrows[i+1][2]:
                double = True
        row_syn.append(double)
    
    # Check if there are resonance structures for reactants
    if row_syn and row_syn[0] == True:
        smiles_dict['reactants_ma'] = smiles_dict.pop('reactants')
        mol_struc["reactants_ma"] = mol_struc.pop('reactants')
        smiles_dict["reactants_mb"] = {"molecules": {}, "arrows": []}
        mol_struc["reactants_mb"] = {}
    
    res = 0
    coordinates = []
    for i in range(len(synthetic_arrows) - 1):
        m = molecule_key(molecules, synthetic_arrows, resonance_arrows, row_syn, i, res, mech_category)
        res = m[-1] +1
        
        if m[0] == 'pass':
            continue
        elif len(m) == 3:
            smiles_dict[m[0]] = {"molecules": {}, "arrows": []}
            mol_struc[m[0]] = {}
            coordinates.append(([m[0]],m[-2]))
        elif len(m) == 4:
            smiles_dict[m[0]] = {"molecules": {}, "arrows": []}
            mol_struc[m[0]] = {}
            smiles_dict[m[1]] = {"molecules": {}, "arrows": []}
            mol_struc[m[1]] = {}
            coordinates.append(([m[0],m[1]],m[-2]))

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
        category = None
        
        # Check if this is a reactant
        if avg_x < synthetic_arrows[0][0]:
            avg_arrow_y = sum([arrow[2] for arrow in synthetic_arrows]) / len(synthetic_arrows) if synthetic_arrows else avg_y
            if row_syn and row_syn[0] == True and avg_y > avg_arrow_y + 50:  # If molecule is below main pathway
                category = "reactants_mb"
            elif row_syn and row_syn[0] == True:
                category = "reactants_ma"
            else:
                category = 'reactants'
            
        # Check if this is a product
        elif not synthetic_arrows or avg_x > synthetic_arrows[-1][1]:
            if not resonance_arrows:
                category = "products"
            elif avg_y < resonance_arrows[0][1]:
                category = "products"

        # Check intermediates
        for x in coordinates:
            x_min =  x[1][0]
            x_max =  x[1][1]
            y =      x[1][2]
            
            if len(x[0]) > 1:
                if x_min < avg_x < x_max and avg_y < y + 50:
                    category = x[0][0]
                elif x_min < avg_x < x_max and avg_y > y + 50:
                    category = x[0][1]
            elif len(x[0]) == 1 and x_min < avg_x < x_max:
                    category = x[0][0]

        
        # If category is found, add molecule to that category
        if category and category in smiles_dict:
            smiles_dict[category]["molecules"][mol_number] = smiles
            molecule_categories[mol_data["a"][0]["i"]] = category
            mol_struc[category][mol_number] = sanitization


    # Process curved arrows
    for arrow in curved_arrows:
        start = arrow['o1']
        end = arrow['o2']
        start_a = None
        end_a = None
        mol_start = None
        mol_end = None
        s = None
        e = None
        
        # Find start atom/bond
        for x in atom_ids:
            if start in atom_ids[x]:
                mol_start = x
                start_a = atom_ids[x].get(start)
                len_start = len(atom_ids[x])
                s = 'atom'
        for x in bond_ids:
            if start in bond_ids[x]:
                mol_start = x
                start_a = bond_ids[x].get(start)
                len_start = len(atom_ids[x])
                s = 'bond'
        
        # Find end atom/bond
        for x in atom_ids:
            if end in atom_ids[x]:
                mol_end = x
                end_a = atom_ids[x].get(end)
                len_end = len(atom_ids[x])
                e = 'atom'
        for x in bond_ids:
            if end in bond_ids[x]:
                mol_end = x
                end_a = bond_ids[x].get(end)
                len_end = len(atom_ids[x])
                e = 'bond'
           
        # Skip if we couldn't find start or end
        if mol_start is None or mol_end is None or s is None or e is None:
            continue
        
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

        # Find which category the start molecule belongs to
        set_key = None
        for key in smiles_dict:
            if mol_start in smiles_dict[key]['molecules'].keys():
                set_key = key
                break
                
        if set_key is not None and start is not None and end is not None:
            smiles_dict[set_key]["arrows"].append({
                "start_molecule": mol_start,
                "start_arrow": start_arrow,
                "substructure_start": SubMol['start_sub'], 
                "end_molecule": mol_end,
                "end_arrow": end_arrow,
                "substructure_end": SubMol['end_sub'],
                "electrons": arrow.get("e")
            })
            
    # Remove empty categories
    keys_to_remove = []
    for key in smiles_dict:
        if not smiles_dict[key]["molecules"] and not smiles_dict[key]["arrows"]:
            keys_to_remove.append(key)
            if key in mol_struc:
                del mol_struc[key]
    
    for key in keys_to_remove:
        del smiles_dict[key]

    if mech_category == 'student':

        return smiles_dict, mol_struc, reaction_arrows

    elif mech_category == 'model':
        return smiles_dict, reaction_arrows



