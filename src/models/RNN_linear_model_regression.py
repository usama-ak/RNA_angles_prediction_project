import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import LSTM, Dense, TimeDistributed



def create_rnn_model(input_shape):
    model = Sequential()
    # LSTM layer with 64 units and return_sequences=True to maintain sequence length
    model.add(LSTM(128, input_shape=input_shape, return_sequences=True))
    model.add(LSTM(128, return_sequences=True))
    model.add(LSTM(64, return_sequences=True))
    # TimeDistributed layer to apply Dense layer to each time step independently
    model.add(TimeDistributed(Dense(1, activation='linear')))  # Outputting an angle per nucleotide
    
    # Compile the model
    model.compile(optimizer='adam', loss='mean_absolute_error', metrics=['mae'])
    
    return model

