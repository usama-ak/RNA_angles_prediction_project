import os
import json
import argparse
import torch
from src.utils.input_preprocess import get_sequences_from_fasta, one_hot_encode_sequence
from src.models.classification.multiple_cl_model import AnglePredictionRNNmulti
from src.models.classification.binary_cl_model import BinaryClassificationRNN

## Setting the paths for the binary and the multi classes pre-trained models
current_dir = os.path.dirname(os.path.abspath(__file__))
multi_model_path = os.path.join(current_dir,'src' , 'models' ,'classification', 'model_multi.pt')
binary_model_path = os.path.join(current_dir,'src' , 'models' ,'classification', 'model_binary.pt')

## Here we define the class AngleHelper that handles angles prediction
class AngleHelper:
    def __init__(self, binary=False):
        self.binary = binary

    def predict(self, in_path: str, out_path: str):
        ## We define the path of the FASTA file then we extract the corresponding sequences
        fasta_file_path = os.path.join(current_dir, in_path)
        sequence = get_sequences_from_fasta(fasta_file_path)
        ## We perform one hot encodening on the sequence
        encoded_sequence, mask = one_hot_encode_sequence(sequence[0])
        encoded_sequence = encoded_sequence.unsqueeze(0)
        mask = mask.unsqueeze(0)
        mask = mask.unsqueeze(0)
        predictions = {}

        if self.binary:
            ## we perform binary angle prediction
            print("Performing binary angle prediction...")
            model = BinaryClassificationRNN()
            model.load_state_dict(torch.load(binary_model_path))
            outputs = model(encoded_sequence, mask)
            _, predicted = torch.max(outputs, 2)
            non_padded_indices = (mask.squeeze(2) == 1).view(-1)
            predicted_non_padded = predicted.view(-1)[non_padded_indices]
            predictions= {
                "sequence": sequence[0],
                "beta angles classes": predicted_non_padded.tolist()
            }

        else:
            ## we perform multiple class angle prediction
            print("Performing multiple class angle prediction...")
            model = AnglePredictionRNNmulti()
            model.load_state_dict(torch.load(multi_model_path))
            outputs = model(encoded_sequence, mask)
            _, predicted = torch.max(outputs, 2)
            non_padded_indices = (mask.squeeze(2) == 1).view(-1)
            predicted_non_padded = predicted.view(-1)[non_padded_indices]
            predictions= {
                "sequence": sequence[0],
                "beta angles classes": predicted_non_padded.tolist()
            }

        ## Then we write predictions to the output file
        with open(out_path, 'w') as output_file:
            json.dump(predictions, output_file, indent=2)
        print(f"Predictions saved to {out_path}")

            

## Here we have the main that will perform angle predictions on a sequence contained in a FASTA file
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Angle prediction for nucleotide sequences")
    parser.add_argument("--binary", action="store_true", help="Perform binary angle prediction")
    args = parser.parse_args()
    
    in_path = os.path.join("data", "sample", "example.fasta")
    out_path = os.path.join(current_dir, "sample.json")
    binary = args.binary

    angle_helper = AngleHelper(binary)
    angle_helper.predict(in_path, out_path)
