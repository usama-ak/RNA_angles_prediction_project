from src.utils.GetSequencesAndAngles import GetSequencesAndAngles
from src.utils.Sequences_padding_and_masking import padding_onehot
from src.utils.dataprepocess import prepare_data, convert_to_tensors
from src.models.Implementation.regression_model import create_rnn_model
from src.models.Evaluation.evaluate_models import print_evaluation_metrics, plot_angles, plot_learning_curves
import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping

train_path = r"C:/Users/usama/Desktop/m2_geniomhe_rna_project-main/data/AngleFilesOutput//"
trainseq_encoded, trainangles_encoded, trainmasks = prepare_data(train_path)

test_path = r"C:/Users/usama/Desktop/m2_geniomhe_rna_project-main/data/AngleFilesTestOutput//"
testseq_encoded, testangles_encoded, testmasks = prepare_data(test_path)

train_seq_tensor, train_angles_tensor, train_masks_tensor = convert_to_tensors(trainseq_encoded, trainangles_encoded, trainmasks)
test_seq_tensor, test_angles_tensor, test_masks_tensor = convert_to_tensors(testseq_encoded, testangles_encoded, testmasks)
'''
print(f"Shape of train_data_X: {train_seq_tensor.shape}")
print(f"Shape of train_data_y: {train_angles_tensor.shape}")
print(f"Shape of train_mask: {train_masks_tensor.shape}")

print(f"Shape of test_data_X: {test_seq_tensor.shape}")
print(f"Shape of test_data_y: {test_angles_tensor.shape}")
print(f"Shape of test_mask: {test_masks_tensor.shape}")'''

early_stopping = EarlyStopping(monitor='val_loss',  # Monitor validation loss
                               patience=5,           # Number of epochs with no improvement after which training will be stopped
                               restore_best_weights=True,  # Restore model weights from the epoch with the best value of the monitored quantity
                               verbose=1)           # Print messages about early stopping to the console


input_shape = (395,4)
rnn_model = create_rnn_model(input_shape)
history = rnn_model.fit(train_seq_tensor,
                        train_angles_tensor,
                        sample_weight=train_masks_tensor,
                        validation_data=(test_seq_tensor, test_angles_tensor, test_masks_tensor),
                        epochs=50,
                        batch_size=32,
                        callbacks=[early_stopping])


'''
evaluation = rnn_model.evaluate(test_seq_tensor, test_angles_tensor, sample_weight=test_masks_tensor)

y_pred = rnn_model.predict(test_seq_tensor)

y_true_np = test_angles_tensor.numpy().flatten()
y_pred_np = y_pred.flatten()

print_evaluation_metrics(evaluation, test_angles_tensor, y_pred)

plot_angles(y_true_np, y_pred_np, title="Expected vs Predicted Angles")

# Assuming you want to test on the first sequence in the test dataset
single_test_sequence_tensor = tf.expand_dims(test_seq_tensor[0], axis=0)

# Predict angles for the single sequence
predicted_angles = rnn_model.predict(single_test_sequence_tensor)

# Print the predicted and real angles
print("Predicted Angles:", predicted_angles.flatten())
print("Real Angles:", test_angles_tensor[0].numpy().flatten())
'''
