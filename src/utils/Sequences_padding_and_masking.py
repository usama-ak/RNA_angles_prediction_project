import tensorflow as tf

def padding_onehot(sequences, angles, max_length):
    # Pad sequences to the maximum length
    padded_sequences = tf.keras.preprocessing.sequence.pad_sequences(sequences, maxlen=max_length, padding='post')
    
    # Calculate the vocabulary size
    vocab_size = 4
    
    # Perform one-hot encoding
    encoded_sequences = tf.one_hot(padded_sequences, depth=vocab_size)

    # Create a mask for the padded values
    mask = tf.cast(tf.math.not_equal(padded_sequences, 0), dtype=tf.float32)

    padded_angles = tf.pad(angles, [[0, max_length - tf.shape(angles)[0]]], constant_values=-1.0)  # Assuming -1 as the padding value
    
    return encoded_sequences, mask, padded_angles
