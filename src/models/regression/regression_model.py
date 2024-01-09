from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, TimeDistributed, GRU, Dropout
from tensorflow.keras.optimizers import Adam

def create_rnn_model(input_shape):
    model = Sequential()
    # LSTM layer with 64 units and return_sequences=True to maintain sequence length
    model.add(GRU(64, input_shape=input_shape, return_sequences=True))
    model.add(Dropout(0.1))
    model.add(LSTM(32, return_sequences=True))
    # TimeDistributed layer to apply Dense layer to each time step independently
    model.add(TimeDistributed(Dense(1, activation='linear')))  # Outputting an angle per nucleotide
    
    custom_learning_rate = 0.001

    optimizer = Adam(learning_rate=custom_learning_rate)

    model.compile(optimizer=optimizer, loss='mean_absolute_error', metrics=['mae'], weighted_metrics=['mae'])
    
    return model 
