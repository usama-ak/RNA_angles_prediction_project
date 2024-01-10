import torch
import torch.nn as nn

class BinaryClassificationRNN(nn.Module):
    def __init__(self, input_size=4, num_layers=4, hidden_size=64, output_size=2, dropout=0.5):
        super(BinaryClassificationRNN, self).__init__()
        self.hidden_size = int(hidden_size.item()) if isinstance(hidden_size, torch.Tensor) else int(hidden_size)
        self.num_layers = int(num_layers.item()) if isinstance(num_layers, torch.Tensor) else int(num_layers)

        # LSTM layers
        self.lstm = nn.LSTM(input_size, self.hidden_size, num_layers=self.num_layers, batch_first=True)
        # Dropout layer
        self.dropout = nn.Dropout(dropout)
        # Fully connected layer
        self.fc = nn.Linear(self.hidden_size, output_size)

    def forward(self, x, mask):
        # Initialize hidden and cell states for LSTM
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        # LSTM forward pass
        output, _ = self.lstm(x, (h0, c0))
        # Apply mask to ignore padded values
        mask_tr = mask.transpose(1, 2) # Squeeze and transpose to move the singleton dimension to the third position
        expanded_mask = mask_tr.expand(-1, -1, output.size(2)) # Expand to match the output tensor
        output = output * expanded_mask # Apply the mask directly
        # Apply dropout
        output = self.dropout(output)
        # Fully connected layer
        output = self.fc(output)
        # output = F.softmax(output, dim=2)
        return output