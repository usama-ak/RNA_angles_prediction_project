import os
import json

## We create a function that will open a json input file
def process_file(input_path, output_path):
    with open(input_path, 'r') as file:
        data = json.load(file)
        ## Here we extract the protein name, the corresponding sequence and the beta angles
        protein_name = list(data.keys())[0]
        sequence = data[protein_name]["sequence"]
        beta_angles = data[protein_name]["angles"]["beta"]
        
        ## Then a dictionnary will be created and the previously extracted information will be added
        output_data = {
            "Protein Name": protein_name,
            "Sequence": sequence,
            "Beta Angles": beta_angles
        }

        ## Writting the result on a json file
        with open(output_path, 'w') as output_file:
            json.dump(output_data, output_file, indent=2)
