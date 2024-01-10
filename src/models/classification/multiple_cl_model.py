import torch
import torch.nn as nn

## We create a class AnglePredictionRNNmulti where we create the multi class instances model that will be used to predict angles
class AnglePredictionRNNmulti(nn.Module):
    def __init__(self, input_size=4, hidden_size=64, output_size=40, num_layers=4, dropout=0.5):
        super(AnglePredictionRNNmulti, self).__init__()
        ## We define hidden_size and number of layers that will be associated to the model
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        # Normalization layer
        #self.norm = nn.LayerNorm(input_size)        
        ## LSTM layer of our model
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers=num_layers, batch_first=True)
        ## Dropout layer of our model
        self.dropout = nn.Dropout(dropout)
        ## Fully connected layer of our model
        self.fc = nn.Linear(hidden_size, output_size)

    def forward(self, x, mask):
        ## Here we initialize the hidden and cell states for LSTM
        h0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, x.size(0), self.hidden_size).to(x.device)
        output, _ = self.lstm(x, (h0, c0))
        ## Then we apply masks to ignore padded values by tansposing and expanding them so that they fit the output tensor
        mask_tr = mask.transpose(1, 2)  
        expanded_mask = mask_tr.expand(-1,-1,output.size(2))  
        output = output * expanded_mask 
        ## Applying dropout
        output = self.dropout(output)
        ## Applying Fully connected layer
        output = self.fc(output)
        return output
