from utils.datapreprocess_cl import prepare_data
import torch
import torch.nn as nn
from models.classification.multiple_cl_model import AnglePredictionRNNmulti
import torch.nn.functional as F
from collections import Counter
from torch.nn.utils import clip_grad_norm_
import os

def categorize_angles(angle):
    if -200 <= angle <= 200:
        return int((angle + 200) / 10)  # 10 degrees per class from -200 to 200

def multi_train(model, train_path, epochs=1000):
    X_train, Y_train, train_masks = prepare_data(train_path)
    Y_train_labels = [ [categorize_angles(angle) for angle in angles] for angles in Y_train ]
    Y_train_labels = torch.tensor(Y_train_labels)
    X_train = X_train.view(-1, X_train.shape[2], X_train.shape[3])
    output_size = 40

    #   Class weight calculation
    Y_train_flat = [item for sublist in Y_train_labels for item in sublist]
    Y_train_flat = [label.item() for label in Y_train_flat]
    # Calculate class frequencies using Counter
    class_counts = Counter(Y_train_flat)
    # Calculate class weights based on the inverse of class frequencies
    total_samples = len(Y_train_flat)
    class_weights = {}
    for i in range(output_size):
        if class_counts[i] != 0:
            class_weights[i] = total_samples / (len(class_counts) * class_counts[i])
        else:
            class_weights[i] = 0  # Assign a weight of 0 to classes with zero occurrences
    class_weights_tensor = torch.tensor([class_weights[i] for i in range(output_size)], dtype=torch.float)

    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

    # Training loop
    for epoch in range(epochs):
        # Set the model to training mode and reset optimizer gradients
        model.train()  
        optimizer.zero_grad()
        
        # Forward pass
        outputs = model(X_train, train_masks)

        # Mask outputs and labels to consider only non-padded values
        transposed_masks = train_masks.transpose(1, 2)  
        expanded_mask = transposed_masks.expand_as(outputs).bool()
        masked_outputs = torch.where(expanded_mask, outputs, torch.zeros_like(outputs))
        masked_labels = torch.where(expanded_mask[..., 0], Y_train_labels, torch.zeros_like(Y_train_labels))

        # Calculate the loss only on non-padded values
        loss = criterion(masked_outputs.view(-1, output_size), masked_labels.long().view(-1))

        # Get predicted classes
        _, predicted = torch.max(outputs, 2)

        # Accuracy calculation using non-padded indices
        non_padded_indices = (train_masks.squeeze(2) == 1).view(-1)
        predicted_non_padded = predicted.view(-1)[non_padded_indices]
        labels_non_padded = Y_train_labels.view(-1)[non_padded_indices]

        correct = (predicted_non_padded == labels_non_padded).sum().item()
        total = non_padded_indices.sum().item()
        accuracy = correct / total

        # Backpropagation and optimization
        loss.backward()
        clip_grad_norm_(model.parameters(), max_norm=1)
        optimizer.step()

        # Print epoch information
        print(f"Epoch: [{epoch+1}/{epochs}] | Loss: {loss:.5f} | Accuracy: {accuracy:.5f}")
    print("Training finished.")


if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, '..', 'data')
    train_path = os.path.join(data_folder, 'AngleFilesOutput')
    model = AnglePredictionRNNmulti()
    multi_train(model, train_path= train_path)









