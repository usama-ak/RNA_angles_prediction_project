from utils.datapreprocess_cl import prepare_data
import torch
import torch.nn as nn
from models.classification.multiple_cl_model import AnglePredictionRNNmulti
import torch.nn.functional as F
from collections import Counter
from torch.nn.utils import clip_grad_norm_
import os

## Here we categorize angles into several classes (40 classes, each representing a 10-degree interval between the range of -200 and 200)
def categorize_angles(angle):
    if -200 <= angle <= 200:
        return int((angle + 200) / 10)  

## Then we use the multiple classes classification model we created before on the train data, by setting a number of epochs equal to 1000
def multi_train(model, train_path, epochs=1000):
    X_train, Y_train, train_masks = prepare_data(train_path)
    ## Here we perform labels categorization
    Y_train_labels = [ [categorize_angles(angle) for angle in angles] for angles in Y_train ]
    Y_train_labels = torch.tensor(Y_train_labels)

    ## Then we reshape the input data
    X_train = X_train.view(-1, X_train.shape[2], X_train.shape[3])
    output_size = 40

    ## Now we calculate the class weights
    Y_train_flat = [item for sublist in Y_train_labels for item in sublist]
    Y_train_flat = [label.item() for label in Y_train_flat]
    ## Then we calculate class frequencies using Counter
    class_counts = Counter(Y_train_flat)
    total_samples = len(Y_train_flat)
    class_weights = {}
    for i in range(output_size):
        if class_counts[i] != 0:
            class_weights[i] = total_samples / (len(class_counts) * class_counts[i]) ## Class weight here depends on the number of occurence of each class
        else:
            class_weights[i] = 0  ## Classes with 0 occurences will have a weight equal to 0
    class_weights_tensor = torch.tensor([class_weights[i] for i in range(output_size)], dtype=torch.float)

    ## Defining the loss function and optimizer
    criterion = nn.CrossEntropyLoss(weight=class_weights_tensor)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

    for epoch in range(epochs):
        model.train()  
        optimizer.zero_grad()
        
        outputs = model(X_train, train_masks)

        ## Then we use torch to define masked outputs and labels so that we take into account only the non-padded values
        transposed_masks = train_masks.transpose(1, 2)  
        expanded_mask = transposed_masks.expand_as(outputs).bool()
        masked_outputs = torch.where(expanded_mask, outputs, torch.zeros_like(outputs))
        masked_labels = torch.where(expanded_mask[..., 0], Y_train_labels, torch.zeros_like(Y_train_labels))

        ## We calculate the loss on the nonpadded values
        loss = criterion(masked_outputs.view(-1, output_size), masked_labels.long().view(-1))

        ## Extracting the predicted classes 
        _, predicted = torch.max(outputs, 2)

        ## Then we use non padded indices to calculate accuracy
        non_padded_indices = (train_masks.squeeze(2) == 1).view(-1)
        predicted_non_padded = predicted.view(-1)[non_padded_indices]
        labels_non_padded = Y_train_labels.view(-1)[non_padded_indices]

        correct = (predicted_non_padded == labels_non_padded).sum().item()
        total = non_padded_indices.sum().item()
        accuracy = correct / total

        ## and finally we perform backpropagation, optimization on the model
        loss.backward()
        clip_grad_norm_(model.parameters(), max_norm=1)
        optimizer.step()

        print(f"Epoch: [{epoch+1}/{epochs}] | Loss: {loss:.5f} | Accuracy: {accuracy:.5f}")
    print("Training finished.")

## Main section of the code that calls the previous function 
if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, '..', 'data')
    train_path = os.path.join(data_folder, 'AngleFilesOutput')
    model = AnglePredictionRNNmulti()
    multi_train(model, train_path= train_path)









