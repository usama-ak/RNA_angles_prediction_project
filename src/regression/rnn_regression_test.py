import tensorflow as tf
from tensorflow.keras.callbacks import EarlyStopping
import os
from utils.evaluate_models import print_evaluation_metrics, plot_angles, plot_learning_curves
from utils.datapreprocess_reg import prepare_data, convert_to_tensors
from models.regression.regression_model import create_rnn_model

current_dir = os.path.dirname(os.path.abspath(__file__))
data_folder = os.path.join(current_dir, '..' ,'..', 'data')

train_path = os.path.join(data_folder, 'AngleFilesOutput')
test_path = os.path.join(data_folder, 'AngleFilesTestOutput')

trainseq_encoded, trainangles_encoded, trainmasks = prepare_data(train_path)
testseq_encoded, testangles_encoded, testmasks = prepare_data(test_path)


train_seq_tensor, train_angles_tensor, train_masks_tensor = convert_to_tensors(trainseq_encoded, trainangles_encoded, trainmasks)
test_seq_tensor, test_angles_tensor, test_masks_tensor = convert_to_tensors(testseq_encoded, testangles_encoded, testmasks)


early_stopping = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True, verbose=1)

input_shape = (395,4)
rnn_model = create_rnn_model(input_shape)
history = rnn_model.fit(train_seq_tensor,
                        train_angles_tensor,
                        sample_weight=train_masks_tensor,
                        validation_data=(test_seq_tensor, test_angles_tensor, test_masks_tensor),
                        epochs=50,
                        batch_size=32,
                        callbacks=[early_stopping])

evaluation = rnn_model.evaluate(test_seq_tensor, test_angles_tensor, sample_weight=test_masks_tensor)
y_pred = rnn_model.predict(test_seq_tensor)
y_true_np = test_angles_tensor.numpy().flatten()
y_pred_np = y_pred.flatten()
print_evaluation_metrics(evaluation, test_angles_tensor, y_pred)
plot_angles(y_true_np, y_pred_np, title="Expected vs Predicted Angles")