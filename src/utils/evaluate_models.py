from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

## The first function here is used to print evaluation metrics (MSE, MAE, Rsquared) when the regression model is used on the train data
def print_evaluation_metrics(evaluation, true_values, predicted_values):
    mse = evaluation[0]
    mae = evaluation[1]
    
    y_true = true_values.numpy().flatten()
    y_pred = predicted_values.flatten()

    ## Calculation of the R squared value
    r_squared = r2_score(y_true, y_pred)

    print("Evaluation MSE:", mse)
    print("Evaluation MAE:", mae)
    print("R-squared:", r_squared)

## This function will plot the learning curves for loss and MAE metrics, respectively.
def plot_learning_curves(history):
    train_loss = history.history['loss']
    val_loss = history.history['val_loss']
    train_mae = history.history['mae']
    val_mae = history.history['val_mae']

    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    plt.plot(train_loss, label='Training Loss')
    plt.plot(val_loss, label='Validation Loss')
    plt.title('Training and Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(train_mae, label='Training MAE')
    plt.plot(val_mae, label='Validation MAE')
    plt.title('Training and Validation MAE')
    plt.xlabel('Epoch')
    plt.ylabel('MAE')
    plt.legend()

    plt.tight_layout()
    plt.show()

## Then this function will print the scatter plot of the expected vs predicted angles values
def plot_angles(y_true, y_pred, title="Expected vs Predicted Angles"):
    plt.figure(figsize=(10, 6))
    
    plt.scatter(y_true, y_pred, alpha=0.5)
    
    ## Here we plot a diagonal line for perfect predictions
    plt.plot([min(y_true), max(y_true)], [min(y_true), max(y_true)], '--', color='red')
    
    plt.title(title)
    plt.xlabel("Expected Angles")
    plt.ylabel("Predicted Angles")
    plt.grid(True)
    plt.show()

