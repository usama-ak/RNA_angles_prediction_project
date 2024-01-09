import os
import json
import numpy as np
from YourModelClass import YourModelClass  
from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot 

class AngleHelper:
    def __init__(self, *args, **kwargs):
        # Initialize any parameters or configurations needed
        self.model_instance = YourModelClass()

    def predict(self, in_path: str, out_path: str):
        # Read the fasta file
        with open(in_path, 'r') as fasta_file:
            fasta_data = fasta_file.readlines()

        sequences = {}
        current_sequence = ""
        for line in fasta_data:
            if line.startswith(">"):
                current_sequence = line[1:].strip()
                sequences[current_sequence] = {"sequence": "", "angles": {"beta": []}}
            else:
                sequences[current_sequence]["sequence"] += line.strip()

        # One-hot encode the sequences using the existing function
        encoded_sequences = {}
        for seq_id, seq_data in sequences.items():
            onehot_encoded = GetSequencesAndAngles(seq_data["sequence"])  
            encoded_sequences[seq_id] = {"sequence": onehot_encoded, "angles": {"beta": []}}

        # Pad and mask the sequences using the existing function
        padded_sequences = padding_onehot([seq_data["sequence"] for seq_data in encoded_sequences.values()])  # Call your existing padding_masking function

        # Use the model for prediction
        for i, (seq_id, seq_data) in enumerate(encoded_sequences.items()):
            predictions = self.model_instance.predict_angles(np.array([padded_sequences[i]]))
            encoded_sequences[seq_id]["angles"]["alpha"] = predictions.tolist()[0]

        # Save the results to a JSON file
        with open(out_path, 'w') as json_file:
            json.dump(encoded_sequences, json_file, indent=2)

if __name__ == "__main__":
    in_path = os.path.join("path_to_your_file", "exemple.fasta")  # Replace with the actual path to your file
    out_path = "sample.json"
    angle_helper = AngleHelper()
    angle_helper.predict(in_path, out_path)

