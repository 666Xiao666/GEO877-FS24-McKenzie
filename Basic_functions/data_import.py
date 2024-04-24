def read_csv(file_path):
    """
    Read a CSV file and return its content as a list of dictionaries.

    Parameters:
    - file_path (str): The path to the CSV file.

    Returns:
    - list of dict: A list where each element is a dictionary representing a row in the CSV file.
    """
    data = []
    with open(file_path, 'r') as file:
        # Read the header
        header = file.readline().strip().split(',')
        for line in file:
            # Split the line into values
            values = line.strip().split(',')
            # Create a dictionary mapping column names to values
            row = {header[i]: values[i] for i in range(len(header))}
            # Append the dictionary to the data list
            data.append(row)
    return data

