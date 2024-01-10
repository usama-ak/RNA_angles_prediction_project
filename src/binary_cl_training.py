from utils.datapreprocess_cl import prepare_data
import torch
import torch.nn as nn
from models.classification.binary_cl_model import BinaryClassificationRNN
import torch.nn.functional as F
import os

## Here we categorize angles into two classes (0 and 1) to perform binary classification 
## Class 0 corresponds to the angles having a value less than 0,  and class 1 the opposite
def categorize_angles(angle):
    if angle <= 0:
        return 0
    else:
        return 1
        
## Then we use the binary classification model we created before on the train data, by setting a number of epochs equal to 1000
def BC_train(model, train_path, epochs=1000):
    X_train, Y_train, train_masks = prepare_data(train_path)
    ## We convert angles to binary labels (0 or 1)
    Y_train_labels = [[categorize_angles(angle) for angle in angles] for angles in Y_train]
    Y_train_labels = torch.tensor(Y_train_labels)
    ## Then reshaping the input data (X_train)
    X_train = X_train.view(-1, X_train.shape[2], X_train.shape[3])
    ## Labels are then flattened (Y_train)
    Y_train_flat = [item for sublist in Y_train_labels for item in sublist]
    ## Then we transform the flattened labels from tensor to a list
    Y_train_flat = [label.item() for label in Y_train_flat]
    ## Now we define the output size,class weights, loss function and the model optimizer
    output_size = 2
    weight_for_class_0 = 3.1 
    class_weights_tensor = torch.tensor([weight_for_class_0] + [1.0] * (output_size - 1), dtype=torch.float)
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()

        outputs = model(X_train, train_masks)

        ## Here, we apply the softmax activation function to the model outputs

        softmax_outputs = F.softmax(outputs, dim=-1)

        ## Then we use torch to defined masked outputs and labels so that we take into account only the non-padded values
        transposed_masks = train_masks.transpose(1, 2)
        expanded_mask = transposed_masks.expand_as(outputs).bool()
        masked_outputs = torch.where(expanded_mask, softmax_outputs, torch.zeros_like(outputs))
        masked_labels = torch.where(expanded_mask[..., 0], Y_train_labels, torch.zeros_like(Y_train_labels))

        ## We calculate the loss on the nonpadded values and perform backpropagation, optimization on the model
        loss = criterion(masked_outputs.view(-1, output_size), masked_labels.long().view(-1))

        accuracy = (torch.argmax(masked_outputs,
                                 dim=-1) == masked_labels.long()).float().sum().item() / masked_labels.numel()

        print(f"Epoch: [{epoch + 1}/{epochs}] | Loss: {loss:.5f} | Accuracy: {accuracy:.5f}")

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
        optimizer.step()

    print("Training finished.")
    
if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, '..', 'data')
    train_path = os.path.join(data_folder, 'AngleFilesOutput')
    model = BinaryClassificationRNNy()
    BC_train(model, train_path= train_path)


