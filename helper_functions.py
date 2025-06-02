from collections import Counter
from difflib import SequenceMatcher

def pathway_difference(values, threshold=50):
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            if abs(values[i] - values[j]) > threshold:
                return True
    return False

def validate_arrows(vectors):
    """
    Validates vectors based on their categories:
    - synthetic/equilibrium vectors should be roughly horizontal and in the same plane (±30)
    - resonance vectors should be roughly vertical and below the synthetic/equilibrium plane
    
    Args:
        vectors: List of tuples (category, x1, x2, y1, y2)
        
    Returns:
        dict: Dictionary with keys as categories and values as booleans indicating validity
    """
    # Initialize result dict to track validity for each category
    result = {
        "synthetic": True,
        "equilibrium": True,
        "resonance": True
    }
    
    # First pass: Find the "plane" of synthetic/equilibrium vectors (average y-value)
    y_values = []
    for vec in vectors:
        category, _, _, y1, y2 = vec
        if category in ["synthetic", "equilibrium"]:
            y_values.extend([y1, y2])
    
    # Calculate the plane (average y-value)
    if y_values:
        plane_y = sum(y_values) / len(y_values)
    else:
        # Default if no synthetic/equilibrium vectors
        plane_y = 0
    
    # Second pass: Check each vector
    for vec in vectors:
        category, x1, x2, y1, y2 = vec
        
        if category in ["synthetic", "equilibrium"]:
            # Check if horizontal (y1 ≈ y2, small difference between y-coordinates)
            is_horizontal = abs(y1 - y2) < 5  # Allowing small difference
            
            # Check if in the same plane (within ±30 of the average y)
            in_plane = abs(y1 - plane_y) <= 30 and abs(y2 - plane_y) <= 30
            
            # Check if x values are different enough to be considered horizontal
            has_x_difference = abs(x1 - x2) > 20  # Reasonable x-difference for horizontal
            
            # Update result for this category if this vector fails
            if not (is_horizontal and in_plane and has_x_difference):
                result[category] = False
                
        elif category == "resonance":
            # Check if vertical (x1 ≈ x2, small difference between x-coordinates)
            is_vertical = abs(x1 - x2) < 50  # Allowing small difference for vertical
            
            # Check if below the plane (y values should be greater than plane_y)
            # Note: In most coordinate systems, "below" means larger y-values
            # If your system is inverted, you may need to reverse this logic
            below_plane = y1 > plane_y + 30 or y2 > plane_y + 30
            
            # Check if y values are different enough to be considered vertical
            has_y_difference = abs(y1 - y2) > 20  # Reasonable y-difference for vertical
            
            # Update result for resonance if this vector fails
            if not (is_vertical and below_plane and has_y_difference):
                result["resonance"] = False
    
    return result

def sort_list_tuple(tuple_list):
    """
    Sorts a list of tuples based on the second value (index 1) in ascending order,
    then extracts only the first value (index 0) from each tuple into a new list.
    
    Args:
        tuple_list: A list of tuples where each tuple has at least 2 values
        
    Returns:
        list: A list containing only the first values from each tuple, sorted by the second value
    """
    # Sort the list of tuples based on the second value (index 1)
    sorted_tuples = sorted(tuple_list, key=lambda x: x[1])
    
    # Extract only the first value (index 0) from each tuple
    result = [item[0] for item in sorted_tuples]
    
    return result

def unique_strings(list_of_lists):
    # Flatten the list and count occurrences
    all_strings = [item for sublist in list_of_lists for item in set(sublist)]
    counts = Counter(all_strings)
    
    # Keep only elements that appear in exactly one sublist
    unique_list = [item for sublist in list_of_lists for item in sublist if counts[item] == 1]
    
    return unique_list

def common_substring(strings):
    if len(strings[0]) > 6 and len(strings[1]) > 6:
        matcher = SequenceMatcher(None, strings[0], strings[1])
        for match in matcher.get_matching_blocks():
            if match.size >= 5:
                return True
    return False

def flatten(xss):
    return [x for xs in xss for x in xs]


def filter_variations_main_res(lst):
    # Extract base items (items without _ma, _mb, _ra, _rb suffix)
    base_items = [item for item in lst if not (item.endswith(("_ma", "_mb", "_ra", "_rb")))]
    
    # Group items by their type (ma/ra and mb/rb)
    a_group = []  # Will contain all items ending with _ma or _ra
    b_group = []  # Will contain all items ending with _mb or _rb
    
    for item in lst:
        if item.endswith("_ma") or item.endswith("_ra"):
            a_group.append(item)
        elif item.endswith("_mb") or item.endswith("_rb"):
            b_group.append(item)
    
    # Generate just two variations:
    # 1. All "a" suffixes (ma + ra)
    # 2. All "b" suffixes (mb + rb)
    
    variations = []
    
    # Variation 1: All "a" items
    a_variation = base_items.copy() + a_group
    a_variation = sorted(a_variation, key=custom_reaction_sort)
    variations.append(a_variation)
    
    # Variation 2: All "b" items
    b_variation = base_items.copy() + b_group
    b_variation = sorted(b_variation, key=custom_reaction_sort)
    variations.append(b_variation)

    return variations

def custom_reaction_sort(item):
    # Determine item category
    if item.startswith('reactants'):
        category = 0  # Reactants come first
    elif item.startswith('intermediates'):
        category = 1  # Intermediates in the middle
        # Extract the intermediate number
        parts = item.split('_')
        if len(parts) >= 2:
            try:
                number = int(parts[1])
                return (category, number)
            except ValueError:
                return (category, 999)
        return (category, 999)
    elif item.startswith('products'):
        category = 2  # Products come last
    else:
        category = 3  # Other items at the end
    
    return (category, 0)

def filter_variations_res(lst):
    base_items = [item for item in lst if not item.endswith(("_a", "_b"))]
    grouped_pairs = {}

    # Group _a and _b pairs by number
    for item in lst:
        if item.endswith("_a") or item.endswith("_b"):
            key = item.rsplit("_", 1)[0]  # Extract base identifier (e.g., "intermediates_1")
            grouped_pairs.setdefault(key, []).append(item)

    # Ensure each key has exactly two elements (_a and _b)
    grouped_pairs = {k: v for k, v in grouped_pairs.items() if len(v) == 2}
    variations = []

    # Generate variations by keeping one element from each pair
    keys = list(grouped_pairs.keys())
    for i in range(2 ** len(keys)):  # Iterate over all possible 2^n combinations
        variation = base_items[:]
        for j, key in enumerate(keys):
            variation.append(grouped_pairs[key][(i >> j) & 1])  # Choose _a or _b based on binary representation
        
        sorted_data = sorted(variation, key=custom_reaction_sort)
        variations.append(sorted_data)
    return variations


def count_trues(d):
    """Count all True values in the first dictionary of a comparison object"""
    true_count = 0
    # Access the first dictionary in the tuple (index 0)
    first_dict = d[0]

    # Iterate through all entries in this dictionary
    for key, value_list in first_dict.items():
        # Each value is a list of tuples: ([bool values], count)
        for bool_tuple in value_list:
            # Count True values in the boolean list
            bool_list = bool_tuple[0]
            true_count += bool_list.count(True)
    
    return true_count
    #return sum(value[0].count(True) for value in d.values())


def list_to_string(lst):
    if len(lst) > 1:
        return ', '.join(map(str, lst[:-1])) + ' and ' + str(lst[-1])
    return ''.join(map(str, lst))

def _split_reaction(reaction_string):
    """Split a reaction string into reactants and products."""
    reactants_str, products_str = reaction_string.split('>>')
    reactants = reactants_str.split('.')
    products = products_str.split('.')
    return reactants, products

def all_lists_empty(dictionary):
    # Iterate through all values in the dictionary
    for value in dictionary.values():
        # Check if the value is a list and is not empty
        if isinstance(value, list) and value:
            return False
    # If we didn't find any non-empty lists, return True
    return True

