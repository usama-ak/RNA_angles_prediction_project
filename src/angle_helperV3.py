import os
import json
import tensorflow as tf
from YourModelClass import YourModelClass
from src.PrepareData import prepare_data

class AngleHelper:
    def __init__(self, *args, **kwargs):
        # Initialize any parameters or configurations needed
        self.model_instance = YourModelClass()

    def predict(self, in_path: str, out_path: str): ## A modifier
        # Call the datapreprocess script
        seq_encoded, angles_encoded, masks = self.prepare_data(in_path)

        # Convert to tensors
        seq_encoded_tensor, angles_encoded_tensor, masks_tensor = self.convert_to_tensors(seq_encoded, angles_encoded, masks)

        # Use the model for prediction
        predictions = self.model_instance.predict_angles([seq_encoded_tensor, masks_tensor])

        # Save the results to a JSON file
        with open(out_path, 'w') as json_file:
            json.dump({"angles": {"beta": predictions.tolist()}}, json_file, indent=2)


if __name__ == "__main__":
    in_path = os.path.join("data", "sample", "exemple.fasta")  # Replace with the actual path to your file
    out_path = "sample.json"
    angle_helper = AngleHelper()
    angle_helper.predict(in_path, out_path)

