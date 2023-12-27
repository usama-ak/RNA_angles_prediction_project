import tensorflow as tf
import os
import glob
import json
import numpy as np
from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot

directory_path = "/home/ali/Downloads/RNA_project/m2_geniomhe_rna_project/data/AngleFilesOutput/"

res=[]
json_files = glob.glob(os.path.join(directory_path, '*.json'))
for json_file in json_files:
        # Process the file using your function
        seq,angles = GetSequencesAndAngles(json_file)
        encoded_result, mask, padded_angles = padding_onehot([seq],angles, max_length)

        # Append the result matrix to the list
        res.append(encoded_result)

masks = [np.any(sequence ==0, axis =1) for sequence in res]

masks = np.array([~mask.astype(int) for mask in masks])



print("Encoded Sequences:")
print(res)
print(masks)
