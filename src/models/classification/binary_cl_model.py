import torch
import torch.nn as nn

## We create a class BinaryClassificationRNN where we create the model that will be used to predict angles
class BinaryClassificationRNN(nn.Module):
    def __init__(self, input_size=4, num_layers=4, hidden_size=64, output_size=2, dropout=0.5):
        super(BinaryClassificationRNN, self).__init__()
        ## We extract hidden_size and number of layers that will be associated to the model
        ## The hyperparameters will be transformed to integers if they are tensors
        self.hidden_size = int(hidden_size.item()) if isinstance(hidden_size, torch.Tensor) else int(hidden_size)
        self.num_layers = int(num_layers.item()) if isinstance(num_layers, torch.Tensor) else int(num_layers)

        ## LSTM layer of our model
        self.lstm = nn.LSTM(input_size, self.hidden_size, num_layers=self.num_layers, batch_first=True)
        ## Dropout layer of our model
        self.dropout = nn.Dropout(dropout)
        ## Fully connected layer of our model
        self.fc = nn.Linear(self.hidden_size, output_size)

    def forward(self, x, mask):
        ## Here we initialize the hidden and cell states for LSTM
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        output, _ = self.lstm(x, (h0, c0))
        ## Then we apply masks to ignore padded values by tansposing and expanding them so that they fit the output tensor
        mask_tr = mask.transpose(1, 2)
        expanded_mask = mask_tr.expand(-1, -1, output.size(2))
        output = output * expanded_mask
        ## Applying dropout
        output = self.dropout(output)
        ## Applying Fully connected layer
        output = self.fc(output)
        # output = F.softmax(output, dim=2)
        return output
