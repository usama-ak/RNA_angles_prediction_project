import json
import tensorflow as tf
import os
import glob
from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot


def convert_to_tensors(seq_encoded, angles_encoded, masks):
    seq_encoded_tensor = tf.convert_to_tensor(seq_encoded)
    seq_encoded_tensor = tf.reshape(seq_encoded_tensor, (seq_encoded_tensor.shape[0], seq_encoded_tensor.shape[2], seq_encoded_tensor.shape[3]))
    angles_encoded_tensor = tf.convert_to_tensor(angles_encoded)
    angles_encoded_tensor = tf.expand_dims(angles_encoded_tensor, axis=-1)
    masks_tensor = tf.convert_to_tensor(masks)
    masks_tensor = tf.transpose(masks_tensor, perm=[0, 2, 1])

    return seq_encoded_tensor, angles_encoded_tensor, masks_tensor

def prepare_data(data_path):
    seq_encoded = []
    seq = []
    angles_encoded = []
    angles = []
    masks = []

    json_files = glob.glob(os.path.join(data_path, '*.json'))
    for json_file in json_files:
        seq_data, beta_angles = GetSequencesAndAngles(json_file)
        seq.append(seq_data)
        angles.append(beta_angles)

    max_length = 395

    for s, a in zip(seq, angles):
        encoded_seq, mask, padded_angles = padding_onehot([s], a, max_length=max_length)
        seq_encoded.append(encoded_seq)
        masks.append(mask)
        angles_encoded.append(padded_angles)

    return seq_encoded, angles_encoded, masks
