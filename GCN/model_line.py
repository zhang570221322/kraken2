import torch
from torch import nn
from torch.nn import functional as F


class FC(nn.Module):

    def __init__(self,  input_dim, output_dim):
        super().__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        print('input dim:', input_dim)
        print('output dim:', output_dim)
        self.out = nn.Sequential(
            nn.Linear(self.input_dim, 512), nn.ReLU(), nn.Linear(512, 128), nn.ReLU(), nn.Linear(128, self.output_dim), nn.LogSigmoid())

    def forward(self, inputs):
        x, adj = inputs
        x = x.squeeze()
        x = self.out(x)
        return x
