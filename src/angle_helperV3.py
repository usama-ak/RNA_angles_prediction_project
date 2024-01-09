import os
import json
import tensorflow as tf
from YourModelClass import YourModelClass
from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot

class AngleHelper:
    def __init__(self, *args, **kwargs):
        # Initialize any parameters or configurations needed
        self.model_instance = YourModelClass()
         
    def convert_to_tensors(self, seq_encoded, angles_encoded, masks):
        seq_encoded_tensor = tf.convert_to_tensor(seq_encoded)
        seq_encoded_tensor = tf.reshape(seq_encoded_tensor, (seq_encoded_tensor.shape[0], seq_encoded_tensor.shape[2], seq_encoded_tensor.shape[3]))
        angles_encoded_tensor = tf.convert_to_tensor(angles_encoded)
        angles_encoded_tensor = tf.expand_dims(angles_encoded_tensor, axis=-1)
        masks_tensor = tf.convert_to_tensor(masks)
        masks_tensor = tf.transpose(masks_tensor, perm=[0, 2, 1])

        return seq_encoded_tensor, angles_encoded_tensor, masks_tensor

    def prepare_data(self, in_path):
        seq_encoded = []
        angles_encoded = []
        masks = []

        seq_data, beta_angles = GetSequencesAndAngles(in_path)
        seq_encoded.append(seq_data)
        angles_encoded.append(beta_angles)

        max_length = 395

        for s, a in zip(seq_encoded, angles_encoded):
            encoded_seq, mask, padded_angles = padding_onehot([s], a, max_length=max_length)
            masks.append(mask)

    def predict(self, in_path: str, out_path: str):
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

