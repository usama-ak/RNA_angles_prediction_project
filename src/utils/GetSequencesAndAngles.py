import json

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
