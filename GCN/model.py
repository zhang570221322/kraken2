import torch
from torch import nn
from torch.nn import functional as F
import pdb

class Arg:
    def __init__(self, num_epochs, learning_rate, weight_decay, batch_size) -> None:
        self.num_epochs = num_epochs
        self.lr = self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.batch_size = batch_size


class GraphConvolution(nn.Module):

    def __init__(self, input_dim, output_dim,
                 dropout=0.):
        super(GraphConvolution, self).__init__()
        self.dropout = dropout
        self.weight = nn.Linear(input_dim, output_dim)
        self.out = nn.Sequential(nn.ReLU(),nn.LayerNorm(output_dim))

    def forward(self, inputs):
        x, adj = inputs
        x = self.weight(x)
        x = torch.bmm(adj, x)
        x = self.out(x) 
        x = F.dropout(x,p=self.dropout, training=self.training)
        return x, adj


class GCN(nn.Module):
    
    def __init__(self,  input_dim, output_dim):
        super(GCN, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        print('input dim:', input_dim)
        print('output dim:', output_dim)
        self.layers = nn.ModuleList()
        self.gcn_dim = [self.input_dim, 8, 8,8,8,8]
        self.num_layer = len(self.gcn_dim)-1
        for i in range(self.num_layer):
            self.layers.append(GraphConvolution(
                self.gcn_dim[i], self.gcn_dim[i+1]))
        self.out = nn.Sequential(
            nn.Linear(8,64), nn.ReLU(), nn.Linear(64,256),nn.ReLU(),nn.Linear(256, self.output_dim))
    def forward(self, inputs):
        x, adj = inputs
        for i in range(len(self.layers)):
            x, adj = self.layers[i]((x, adj))
        # max pool
        x= torch.max(x, dim=1)[0]
        # pdb.set_trace()
        x = self.out(x)
        return x

