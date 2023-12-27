import os
import json
from src.utils.Generating_beta_files import process_file

input_directory = "/home/ali/Downloads/RNA_project/m2_geniomhe_rna_project/data/TestSetOutput"
output_directory = "/home/ali/Downloads/RNA_project/m2_geniomhe_rna_project/data/AngleFilesTestOutput"


if not os.path.exists(output_directory):
    os.makedirs(output_directory)

    # List all files in the input directory
    files = [f for f in os.listdir(input_directory) if f.endswith(".pdb.json")]

    # Process each file
    for file_name in files:
        input_path = os.path.join(input_directory, file_name)
        output_path = os.path.join(output_directory, f"{file_name}_output.json")
        process_file(input_path, output_path)
        


