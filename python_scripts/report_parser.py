import re
import pandas as pd

def parse_file(file_path):
    """
    Parses a file with fixed and variable key-value pairs into a dictionary and optionally a DataFrame.

    Args:
        file_path (str): The path to the file to be parsed.

    Returns:
        dict: A dictionary containing the parsed data.
        pd.DataFrame: A DataFrame representation of the parsed data.
    """
    parsed_data = {}

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()
            if ":" in line:
                key, value = line.split(":")
                parsed_data[key.strip()] = float(value.strip()) if '.' in value else int(value.strip())

    # Convert the dictionary to a DataFrame
    df = pd.DataFrame(parsed_data.items(), columns=["Key", "Value"])

    return parsed_data, df
