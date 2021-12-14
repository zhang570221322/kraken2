import torch
from torch import nn
from torch.nn import functional as F


class GraphConvolution(nn.Module):

    def __init__(self, input_dim, output_dim,
                 dropout=0.,
                 bias=True,
                 activation=F.relu):
        super(GraphConvolution, self).__init__()

        self.dropout = dropout
        self.bias = bias
        self.activation = activation
        self.weight = nn.Parameter(torch.randn(
            input_dim,  output_dim))
        if bias:
            self.bias = nn.Parameter(torch.zeros(
                output_dim))

    def forward(self, inputs):
        x, adj = inputs
        if self.training:
            x = F.dropout(x, self.dropout)
        x = torch.matmul(adj.float(), x)
        x = torch.matmul(x, self.weight)
        if self.bias is not None:
            x += self.bias
        x = self.activation(x)
        return x


class GCN(nn.Module):

    def __init__(self,  input_dim, output_dim):
        super(GCN, self).__init__()

        self.input_dim = input_dim
        self.output_dim = output_dim
        print('input dim:', input_dim)
        print('output dim:', output_dim)
        self.layers = nn.Sequential(GraphConvolution(self.input_dim, 16,
                                                     dropout=0.5),
                                    GraphConvolution(16,   64,
                                                     dropout=0.5),
                                    GraphConvolution(64,  256,
                                                     dropout=0.5),
                                    )
        self.out = nn.Sequential(
            nn.Linear(256, 32), nn.ReLU(), nn.Linear(32, 1))

    def forward(self, inputs):
        x, adj = inputs
        x = self.layers((x, adj))
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
