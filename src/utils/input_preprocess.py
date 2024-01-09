from Bio import SeqIO

def get_sequences_from_fasta(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

import torch

def one_hot_encode_sequence(sequence):
    # Define the mapping for one-hot encoding
    mapping = {'A': 1, 'C': 2, 'G': 3, 'U': 4}
    
    # Convert the sequence to a list of integers using the mapping
    sequence_int = [mapping.get(base, -1) for base in sequence]
    # Pad sequence to the maximum length
    max_length = 395
    padded_sequence_int = sequence_int + [0] * (max_length - len(sequence_int))

    # Perform one-hot encoding manually
    vocab_size = 4
    encoded_sequence = torch.zeros(max_length, vocab_size)
    mask = torch.zeros(max_length)

    for i, base_idx in enumerate(padded_sequence_int):
        if base_idx != 0:  # Ignore padding value
            encoded_sequence[i, base_idx - 1] = 1  # Assuming the base index starts from 1
            mask[i] = 1  # Set mask to 1 for real nucleotide

    return encoded_sequence, mask

