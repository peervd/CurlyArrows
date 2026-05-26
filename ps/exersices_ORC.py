import csv
import os
from rdkit import Chem

def csv_path(file_name: str) -> str:
    """Return absolute path to a CSV inside the ps/ directory (next to this file)."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, file_name)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find CSV: {path}")
    return path

def dict_exe(file_name: str):
    """
    Reads a CSV file with ';' as a separator and converts it into a dictionary.
    - The first column stores a dictionary with row numbers as keys and first column values.
    - Other columns store lists with column names as keys.
    """

    # Determine the absolute path to the ps/ folder
    base_dir = os.path.dirname(os.path.abspath(__file__))  
    file_path = os.path.join(base_dir, file_name)

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Could not find CSV: {file_path}")

    data_dict = {}
    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.reader(file, delimiter=";")
        headers = next(reader)  # Get column names

        # Initialize dictionary with empty lists for all columns except the first
        data_dict[headers[0]] = {}  # First column stores a dictionary
        for header in headers[1:]:
            data_dict[header] = []

        # Populate dictionary with data
        for index, row in enumerate(reader, start=1):
            data_dict[headers[0]][index] = row[0]  # First column as dict
            for i in range(1, len(headers)):  
                data_dict[headers[i]].append(row[i])

    return data_dict


def exersices(number):
    exersice_info = dict_exe('exercise_model_answers.csv')
    if number not in exersice_info[list(exersice_info.keys())[0]].keys():
        print('exersice does not exist, check if you have the correct exersice number')
    code = {}
    for i, key in enumerate(exersice_info):
        if i == 0:
            code[key] = exersice_info[key][number]
        else:
            code[key] = exersice_info[key][number-1]

    return code

def get_acid_base(model,acid_base_file):
    allowed_ab = load_acid_base(acid_base_file)
    all_comp = [sublist for x in model.keys() for sublist in list(model[x]['molecules'].values())]
    acid = ['[H+]']
    base = ['[BaH2]','[BaH3-]']
    if 'O'in all_comp:
        acid.append('[OH3+]')
    for a,b in allowed_ab:
        try:
            a_mol = Chem.MolFromSmiles(a)
            a_minH_mol = Chem.RemoveHs(a_mol)       
            a_minH = Chem.MolToSmiles(a_minH_mol)
        except:
            a_minH = 'NA_a'
        if a in all_comp or a_minH in all_comp and a != '[H+]':
            acid.append(a)
            if a_minH != 'NA_a':
                acid.append(a_minH)
        try:
            b_mol = Chem.MolFromSmiles(b)
            b_minH_mol = Chem.RemoveHs(b_mol)       
            b_minH = Chem.MolToSmiles(b_minH_mol)
        except:
            b_minH = 'NA_b'
        
        if b in all_comp:
            base.append(b)
            if b_minH != 'NA_b':
                base.append(b_minH)
    if not base:
        base.append('no_base')
    return acid,base
    
    
def load_acid_base(file_name: str):
    """
    Reads a CSV file from the ps/ folder and converts it into a list of pairs.
    - Skips the header row.
    - Each row becomes [row[0], row[1]].
    """

    # Get absolute path inside ps/
    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, file_name)

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Could not find CSV: {file_path}")

    acid_base = []
    with open(file_path, mode="r", encoding="utf-8") as file:
        reader = csv.reader(file)
        next(reader, None)  # Skip header row safely
        for row in reader:
            if len(row) >= 2:  # Safety: only take valid rows
                acid_base.append([row[0], row[1]])

    return acid_base

def acid_base(file_path):
    """
    Reads a CSV file and converts it into a dictionary.
    The keys are the row numbers (starting from 1), and the values are from the second column.
    """
    file_path = csv_path(file_path)
    
    acid_base = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            a_b = [row[0],row[1]]
            acid_base.append(a_b) 
    
    acid = []
    for a in acid_base:
        acid.append(a[0])
    base = []
    for b in acid_base:
        base.append(b[1])
    return [acid,base]


def get_concept_tags(file_path, reaction_number):

    """
    Load reaction data from a CSV file for a specific reaction number.
    
    Args:
        file_path (str): Path to the CSV file
        reaction_number (int): The reaction number to filter for
    
    Returns:
        dict: Dictionary containing reaction type and steps with their details
    """
    file_path = csv_path(file_path)
    result = {'reaction_type': None,'steps': {}}
    
    with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile,delimiter="/")
        
        for row in reader:
            # Check if this row matches the reaction number
            if int(row['reaction']) == reaction_number:
                # Set reaction type if not already set
                if result['reaction_type'] is None:
                    result['reaction_type'] = row['description']
                    result['reactants'] = row['reactants']
                    result['product'] = row['product']
                # Add step information
                step_num = int(row['step'])
                result['steps'][step_num] = {
                    'category': row['category'],
                    'tags_1': row['concept_1'],
                    'tags_2': row['concept_2']}

    # Return None if no matching reaction was found
    if result['reaction_type'] is None:
        return None
    else:    
        return result