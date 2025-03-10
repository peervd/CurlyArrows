import csv

def dict_exe(file_path):
    """
    Reads a CSV file and converts it into a dictionary.
    The keys are the row numbers (starting from 1), and the values are from the second column.
    """
    data_dict = {}
    
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        for index, row in enumerate(reader, start=1):
            data_dict[index] = row[0]  # Store the second column as value
    
    return data_dict


def exersices(number):
    java_code = dict_exe('exercise_model_answers.csv')
    if number not in java_code.keys():
        print('exersice does not exist, check if you have the correct exersice number')
        
    code = java_code[number]
    return code

