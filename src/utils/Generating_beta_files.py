import os
import json

def process_file(input_path, output_path):
    with open(input_path, 'r') as file:
        data = json.load(file)
        protein_name = list(data.keys())[0]
        sequence = data[protein_name]["sequence"]
        beta_angles = data[protein_name]["angles"]["beta"]
        
        # Create a dictionary with the required information
        output_data = {
            "Protein Name": protein_name,
            "Sequence": sequence,
            "Beta Angles": beta_angles
        }

        # Write the information to a JSON file
        with open(output_path, 'w') as output_file:
            json.dump(output_data, output_file, indent=2)
