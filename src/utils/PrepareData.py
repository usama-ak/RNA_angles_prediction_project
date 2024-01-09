import torch
from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot

def prepare_data(data_path, max_length=395):
    seq_encoded = []
    angles_encoded = []
    masks = []

    json_files = glob.glob(os.path.join(data_path, '*.json'))
    for json_file in json_files:
        seq_data, beta_angles = GetSequencesAndAngles(json_file)
        encoded_seq, mask, padded_angles = padding_onehot([seq_data], beta_angles, max_length=max_length)
        seq_encoded.append(encoded_seq)
        masks.append(mask)
        angles_encoded.append(padded_angles)

    # Convert lists to tensors before returning
    seq_encoded = torch.stack(seq_encoded)
    angles_encoded = torch.stack(angles_encoded)
    masks = torch.stack(masks)

    return seq_encoded, angles_encoded, masks
