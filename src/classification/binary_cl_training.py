from src.preprocessing.dataprepocess_cl import prepare_data
import torch
import torch.nn as nn
from src.models.classification.binary_cl_model import BinaryClassificationRNN
import torch.nn.functional as F

def categorize_angles(angle):
    if angle <= 0:
        return 0
    else:
        return 1

def BC_train(model, train_path="path_to_AngleFilesOutput", epochs=1000):
    X_train, Y_train, train_masks = prepare_data(train_path)
    Y_train_labels = [[categorize_angles(angle) for angle in angles] for angles in Y_train]
    Y_train_labels = torch.tensor(Y_train_labels)
    X_train = X_train.view(-1, X_train.shape[2], X_train.shape[3])
    Y_train_flat = [item for sublist in Y_train_labels for item in sublist]
    # Transform Y_train_flat from a tensor to a list
    Y_train_flat = [label.item() for label in Y_train_flat]
    output_size = 2
    weight_for_class_0 = 3.1 
    class_weights_tensor = torch.tensor([weight_for_class_0] + [1.0] * (output_size - 1), dtype=torch.float)
    # Loss function and optimizer
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        # Forward pass
        outputs = model(X_train, train_masks)

        # Apply softmax activation
        softmax_outputs = F.softmax(outputs, dim=-1)

        # Mask outputs and labels to consider only non-padded values
        transposed_masks = train_masks.transpose(1, 2)
        expanded_mask = transposed_masks.expand_as(outputs).bool()
        masked_outputs = torch.where(expanded_mask, softmax_outputs, torch.zeros_like(outputs))
        masked_labels = torch.where(expanded_mask[..., 0], Y_train_labels, torch.zeros_like(Y_train_labels))

        # Calculate the loss only on non-padded values
        loss = criterion(masked_outputs.view(-1, output_size), masked_labels.long().view(-1))

        # Print epoch information
        accuracy = (torch.argmax(masked_outputs,
                                 dim=-1) == masked_labels.long()).float().sum().item() / masked_labels.numel()

        print(f"Epoch: [{epoch + 1}/{epochs}] | Loss: {loss:.5f} | Accuracy: {accuracy:.5f}")

        # Backpropagation and optimization
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
        optimizer.step()

    print("Training finished.")
if __name__ == "__main__":
    # Initialize the model
    model = AnglePredictionRNN_binary()
    BC_train(model)


