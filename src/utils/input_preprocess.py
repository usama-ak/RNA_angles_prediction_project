from Bio import SeqIO

## First we define a function that will extract the sequences from a FASTA file using SeqIO
def get_sequences_from_fasta(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

import torch

## Then we have a second function for which one-hot enconding, padding and masking will be done for the sequence
def one_hot_encode_sequence(sequence):
    ## Here we define the mapping for one-hot encoding
    mapping = {'A': 1, 'C': 2, 'G': 3, 'U': 4}
    
    ## Then, converting the sequence to a list of integers using the mapping
    sequence_int = [mapping.get(base, -1) for base in sequence]
    ## Padding the sequence to a maximum length so that it fits with the input size of the model
    max_length = 395
    padded_sequence_int = sequence_int + [0] * (max_length - len(sequence_int))

    vocab_size = 4
    encoded_sequence = torch.zeros(max_length, vocab_size)
    mask = torch.zeros(max_length) ## Generating masks for the sequence

    for i, base_idx in enumerate(padded_sequence_int):
        if base_idx != 0:  
            encoded_sequence[i, base_idx - 1] = 1  
            mask[i] = 1  ## Here we set the mask value to 1 for real nucleotides

    return encoded_sequence, mask

