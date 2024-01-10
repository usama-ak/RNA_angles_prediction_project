import torch
import torch.nn as nn

class AnglePredictionRegression(nn.Module):
    def __init__(self, input_size=4, hidden_size=64, output_size=1, num_layers=4, dropout=0.5):
        super(AnglePredictionRegression, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers=num_layers, batch_first=True)
        self.dropout = nn.Dropout(dropout)
        self.fc = nn.Linear(hidden_size, output_size)
        
    def forward(self, x, mask):
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        
        output, _ = self.lstm(x, (h0, c0))

        # Apply mask to ignore padded values
        if mask.shape[-2:] == torch.Size([1, 395]):
            mask_tr = mask.transpose(1, 2)
        else:
            mask_tr = mask
        
        expanded_mask = mask_tr.expand(-1, -1, output.size(2)) 
        output = output * expanded_mask
        
        output = self.dropout(output)
        output = self.fc(output)
        
        return output
