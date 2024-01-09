import torch
import torch.nn as nn

class AnglePredictionRNN(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, num_layers=4, dropout=0.5):
        super(AnglePredictionRNN, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        # Normalization layer
        #self.norm = nn.LayerNorm(input_size)        
        # LSTM layers
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers=num_layers, batch_first=True)
        # Dropout layer
        self.dropout = nn.Dropout(dropout)
        # Fully connected layer
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x, mask):
        # Initialize hidden and cell states for LSTM
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        # LSTM forward pass
        output, _ = self.lstm(x, (h0, c0))
        # Apply mask to ignore padded values
        mask_tr = mask.transpose(1, 2)  
        expanded_mask = mask_tr.expand(-1,-1,output.size(2))  
        output = output * expanded_mask 
        # Apply dropout
        output = self.dropout(output)
        # Fully connected layer
        output = self.fc(output)
        return output