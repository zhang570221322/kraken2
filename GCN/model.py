import torch
from torch import nn
from torch.nn import functional as F


class GraphConvolution(nn.Module):

    def __init__(self, input_dim, output_dim,
                 dropout=0.,
                 activation=F.relu):
        super(GraphConvolution, self).__init__()

        self.dropout = dropout
        self.activation = activation
        self.weight = nn.Linear(input_dim, output_dim)

    def forward(self, inputs):
        x, adj = inputs
        x = self.weight(x)
        x = torch.matmul(adj, x)
        x = F.dropout(self.activation(x))
        return x, adj


class GCN(nn.Module):

    def __init__(self,  input_dim, output_dim):
        super(GCN, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        print('input dim:', input_dim)
        print('output dim:', output_dim)
        self.layers = nn.ModuleList()
        self.num_layer = 3
        self.gcn_dim = [self.input_dim, 16, 32, 64]
        for i in range(self.num_layer):
            self.layers.append(GraphConvolution(
                self.gcn_dim[i], self.gcn_dim[i+1]))
        self.out = nn.Sequential(
            nn.Linear(64, 32), nn.ReLU(), nn.Linear(32, self.output_dim))

    def forward(self, inputs):
        x, adj = inputs
        for i in range(len(self.layers)):
            x, adj = self.layers[i]((x, adj))
        x = self.out(x)
        return x

    def l2_loss(self):
        layer = self.layers.children()
        layer = next(iter(layer))
        loss = None
        for p in layer.parameters():
            if loss is None:
                loss = p.pow(2).sum()
            else:
                loss += p.pow(2).sum()
        return loss
