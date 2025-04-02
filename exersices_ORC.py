import csv

def dict_exe(file_path):
    """
    Reads a CSV file with ';' as a separator and converts it into a dictionary.
    - The first column stores a dictionary with row numbers as keys and first column values.
    - Other columns store lists with column names as keys.
    """
    data_dict = {}

    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter=';')
        headers = next(reader)  # Get column names
        
        # Initialize dictionary with empty lists for all columns except the first
        data_dict[headers[0]] = {}  # First column stores a dictionary
        for header in headers[1:]:
            data_dict[header] = []

        # Populate dictionary with data
        for index, row in enumerate(reader, start=1):
            data_dict[headers[0]][index] = row[0]  # First column stored as a dictionary
            for i in range(1, len(headers)):  # Other columns stored as lists
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


def load_acid_base(file_path):
    """
    Reads a CSV file and converts it into a dictionary.
    The keys are the row numbers (starting from 1), and the values are from the second column.
    """
    acid_base = []
    
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            a_b = [row[0],row[1]]
            acid_base.append(a_b) 
    
    return acid_base