import os
import argparse
from input_preprocess import get_sequences_from_fasta, one_hot_encode_sequence

class AngleHelper:
    def __init__(self, binary=False):
        self.binary = binary

    def predict(self, in_path: str, out_path: str):
        sequences = self.get_sequences_from_fasta(in_path)
        predictions = {}
        if self.binary:
            # Implement logic for binary angle prediction
            print("Performing binary angle prediction...")
        else:
            # Implement logic for multiple class angle prediction
            print("Performing multiple class angle prediction...")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Angle prediction for nucleotide sequences")
    parser.add_argument("input", help="Path to the input FASTA file")
    parser.add_argument("output", help="Path to the output JSON file")
    parser.add_argument("--binary", action="store_true", help="Perform binary angle prediction")

    args = parser.parse_args()
    
    in_path = args.input
    out_path = args.output
    binary = args.binary

    angle_helper = AngleHelper(binary)
    angle_helper.predict(in_path, out_path)
