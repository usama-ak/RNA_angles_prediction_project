from utils.datapreprocess_cl import prepare_data
import torch
import torch.nn as nn
from models.regression.regression_model import AnglePredictionRegression
import os


def regression_train(model, train_path, epochs=30):
    X_train, Y_train, train_masks = prepare_data(train_path)
    X_train = X_train.view(-1, X_train.shape[2], X_train.shape[3])

    criterion = nn.MSELoss()  # Using Mean Squared Error for regression
    optimizer = torch.optim.Adam(model.parameters(), lr=0.1)
    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        outputs = model(X_train, train_masks)
        criterion = nn.MSELoss()
        # Masking to handle padded values
        if train_masks.shape[-2:] == torch.Size([1, 395]):
            train_masks = train_masks.transpose(1, 2)
        else:
            train_masks = train_masks    
        masked_outputs = outputs[train_masks.squeeze(-1).bool()]
        masked_labels = Y_train[train_masks.squeeze(-1).bool()]
        # Calculate loss on non-padded values
        loss = criterion(masked_outputs.view(-1), masked_labels.view(-1))
        # Calculate metrics (e.g., MAE)
        mae = torch.abs(masked_outputs - masked_labels).mean().item()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1)
        optimizer.step()
        print(f"Epoch: [{epoch + 1}/{epochs}] | Loss: {loss.item()} | Mean Absolute Error: {mae}")
    print("Training finished.")


def regression_test(model, test_path):
    X_test, Y_test, test_masks = prepare_data(test_path)
    X_test = X_test.view(-1, X_test.shape[2], X_test.shape[3])

    # Put the model in evaluation mode
    model.eval()

    with torch.no_grad():
        # Get predictions for the test set
        test_outputs = model(X_test, test_masks)
        # Masking to handle padded values
        if test_masks.shape[-2:] == torch.Size([1, 395]):
            test_masks = test_masks.transpose(1, 2)
        else:
            test_masks = test_masks    
        masked_test_outputs = test_outputs[test_masks.squeeze(-1).bool()]
        masked_test_labels = Y_test[test_masks.squeeze(-1).bool()]

        # Calculate MSE and MAE
        mse = nn.MSELoss()(masked_test_outputs.view(-1), masked_test_labels.view(-1)).item()
        mae = torch.abs(masked_test_outputs - masked_test_labels).mean().item()

        # Calculate R-squared
        mean_y = torch.mean(masked_test_labels)
        ss_total = torch.sum((masked_test_labels - mean_y) ** 2)
        ss_residual = torch.sum((masked_test_labels - masked_test_outputs) ** 2)
        r_squared = 1 - (ss_residual / ss_total)

        # Plotting expected vs predicted angles for the test dataset
        import matplotlib.pyplot as plt

        plt.figure(figsize=(8, 6))
        plt.scatter(masked_test_labels, masked_test_outputs, alpha=0.5)
        plt.title('Expected vs Predicted Angles')
        plt.xlabel('Expected Angles')
        plt.ylabel('Predicted Angles')
        plt.grid()
        plt.show()

        # Print MSE and MAE
        print(f"Mean Squared Error (MSE): {mse}")
        print(f"Mean Absolute Error (MAE): {mae}")
        print(f"R-squared (R²): {r_squared}")


if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(current_dir, '..', 'data')
    train_path = os.path.join(data_folder, 'AngleFilesOutput')
    test_path = os.path.join(data_folder, 'AngleFilesTestOutput')
    model = AnglePredictionRegression()
    regression_train(model, train_path=train_path)
    regression_test(model, test_path=test_path)
