import json
import tensorflow as tf
import os
import glob

def GetSequencesAndAngles(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
        sequence = data["Sequence"]
        angles = data["Beta Angles"]

    # Create a mapping for one-hot encoding
    mapping = {'A': 1, 'C': 2, 'G': 3, 'U': 4}

    # Convert the sequence to a list of integers using the mapping
    sequence_int = [mapping.get(base, -1) for base in sequence]

    # Convert beta angles to float
    angles_float = [float(angle) for angle in angles if angle.strip()]

    return sequence_int,angles_float

def padding_onehot(sequences, angles, max_length):
    # Pad sequences to the maximum length
    padded_sequences = tf.keras.preprocessing.sequence.pad_sequences(sequences, maxlen=max_length, padding='post')
    
    # Calculate the vocabulary size
    vocab_size = 4
    
    # Perform one-hot encoding
    encoded_sequences = tf.one_hot(padded_sequences, depth=vocab_size)

    # Create a mask for the padded values
    mask = tf.cast(tf.math.not_equal(padded_sequences, 0), dtype=tf.float32)

    padded_angles = tf.pad(angles, [[0, max_length - tf.shape(angles)[0]]], constant_values=-1.0)  # Assuming -1 as the padding value
    
    return encoded_sequences, mask, padded_angles

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