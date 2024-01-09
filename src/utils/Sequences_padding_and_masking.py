import tensorflow as tf
import torch


def padding_onehot(sequences, angles, max_length):
    # Pad sequences to the maximum length
    padded_sequences = [seq + [0] * (max_length - len(seq)) for seq in sequences]
    padded_sequences_tensor = torch.tensor(padded_sequences)  # Convert to PyTorch tensor
    
    # Perform one-hot encoding manually
    vocab_size = 4
    encoded_sequences = torch.zeros(len(padded_sequences), max_length, vocab_size)
    for i, seq in enumerate(padded_sequences):
        for j, base in enumerate(seq):
            if base != 0:  # Ignore padding value
                encoded_sequences[i, j, base - 1] = 1  # Assuming the base index starts from 1

    # Create a mask for the padded values
    mask = (padded_sequences_tensor != 0).float()

    # Pad angles manually
    padded_angles = torch.tensor(angles + [-1] * (max_length - len(angles)))

    return encoded_sequences, mask, padded_angles
