# RNA Angle Prediction with RNN (Binary/Multiple Class)

This code performs RNA angle prediction using RNN-based models for both binary and multiple class classification.

## Requirements
The requitements are specified in the `requirements.txt` file.

## Usage

### Running Predictions

1. Clone the repository:
    ```bash
    git clone https://github.com/usama-ak/RNA_angles_prediction_project.git
    ```

2. Navigate to the directory:
    ```bash
    cd RNA_angles_prediction_project
    ```

3. Run predictions:
    ```bash
    python angle_helper_vf.py --binary
    ```
   For multiple class prediction, exclude the `--binary` flag.

### Input/Output

- Input: Fasta file with RNA sequence.
- Output: JSON file containing sequence and predicted beta angle classes.

### Output Classes

- **Binary Classification**:
  - `0`: Represents negative angles.
  - `1`: Represents positive angles.

- **Multiple Class Classification**:
  - Total Classes: 40
  - Each class represents a 10-degree range from -200 to 200.

## Files

### `angle_helper_vf.py`

- `AngleHelper`: Python class to handle angle predictions.
- `predict`: Method for predicting angles based on input sequences.
- `--binary`: Flag to perform binary angle prediction.

### Folders

- `src/models/classification/`: Contains RNN models for binary and multiple class classification.
- `data/sample/example.fasta`: Sample input file with RNA sequences.
- `sample.json`: Output file with predicted angles in JSON format.


## Authors

- [Usama AKHTAR](https://github.com/usama-ak)
- [Ali YOUCHA](https://github.com/MrAli1582)
- [Abdelouahab ELKOUADI]()
