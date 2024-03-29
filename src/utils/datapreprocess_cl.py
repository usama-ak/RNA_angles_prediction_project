import json
import torch
import os
import glob

## This first function is dedicated to the extraction of the sequence and the beta angles from a JSON file
def GetSequencesAndAngles(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
        sequence = data["Sequence"]
        angles = data["Beta Angles"]
    ## Here we define the mapping for one-hot encoding
    mapping = {'A': 1, 'C': 2, 'G': 3, 'U': 4}
    ## Then, converting the sequence to a list of integers using the mapping
    sequence_int = [mapping.get(base, -1) for base in sequence]
    ## And we convert beta angles to float
    angles_float = [float(angle) for angle in angles if angle.strip()]
    return sequence_int,angles_float

## The second function is dedicated to performing padding
def padding_onehot(sequences, angles, max_length):
    ## Padding the sequence to a maximum length
    padded_sequences = [seq + [0] * (max_length - len(seq)) for seq in sequences]
    padded_sequences_tensor = torch.tensor(padded_sequences) 
    
    vocab_size = 4
    encoded_sequences = torch.zeros(len(padded_sequences), max_length, vocab_size)
    for i, seq in enumerate(padded_sequences):
        for j, base in enumerate(seq):
            if base != 0:  
                encoded_sequences[i, j, base - 1] = 1  

    ## Generating a mask for the padded values
    mask = (padded_sequences_tensor != 0).float()

    ## Padding angles manually
    padded_angles = torch.tensor(angles + [-1] * (max_length - len(angles)))

    return encoded_sequences, mask, padded_angles

## Then the third function will be dedicated to data preparation for the multiple classes classification approach
def prepare_data(data_path, max_length=395):
    seq_encoded = []
    angles_encoded = []
    masks = []

    json_files = glob.glob(os.path.join(data_path, '*.json'))
    for json_file in json_files:
        seq_data, beta_angles = GetSequencesAndAngles(json_file) ## We extract the sequence and beta angles using GetSequencesAndAngles function
        encoded_seq, mask, padded_angles = padding_onehot([seq_data], beta_angles, max_length=max_length) ## We generate the padded angles and masks using padding_onehot function
        seq_encoded.append(encoded_seq)
        masks.append(mask)
        angles_encoded.append(padded_angles)

    ## Converting lists to tensors
    seq_encoded = torch.stack(seq_encoded)
    angles_encoded = torch.stack(angles_encoded)
    masks = torch.stack(masks)

    return seq_encoded, angles_encoded, masks
