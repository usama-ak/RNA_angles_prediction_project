from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, TimeDistributed, GRU, Dropout
from tensorflow.keras.optimizers import Adam

## We create a regression model 
def create_rnn_model(input_shape):
    model = Sequential()
    ## Then we add three layers to the model : GRU, Dropout (to drop some of the neurons) and LSTM
    model.add(GRU(64, input_shape=input_shape, return_sequences=True))
    model.add(Dropout(0.1))
    model.add(LSTM(32, return_sequences=True))
    ## We also add a TimeDistributed layer to apply Dense layer to each time step independently (so that we can output an angle value per nucleotide)
    model.add(TimeDistributed(Dense(1, activation='linear'))) 

    ## Then we define hyperparameters (learning rate, optimizer)
    custom_learning_rate = 0.001

    optimizer = Adam(learning_rate=custom_learning_rate)

    model.compile(optimizer=optimizer, loss='mean_absolute_error', metrics=['mae'], weighted_metrics=['mae'])
    
    return model 
