import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, TopKPooling
from torch_geometric.nn import global_mean_pool as gap
from torch_geometric.nn import global_max_pool as gmp
import pdb


class Net(nn.Module):
    def __init__(self, dataset, hidden_channels) -> None:
        super().__init__()
        torch.manual_seed(12345)
        hcs = hidden_channels
        self.conv1 = GCNConv(dataset.num_node_features, hcs)
        self.pool1 = TopKPooling(hidden_channels, ratio=0.5)
        self.conv2 = GCNConv(hcs, hcs)
        self.pool2 = TopKPooling(hidden_channels, ratio=0.5)
        self.out = nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.3, training=self.training),
            nn.Linear(hcs*2, hcs*4), nn.ReLU(), nn.Dropout(p=0.3,
                                                           training=self.training),
            nn.Linear(hcs*4, hcs*8), nn.ReLU(), nn.Dropout(p=0.3,
                                                           training=self.training),
            nn.Linear(hcs*8, dataset.num_classes), nn.Sigmoid())

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        x = self.conv1(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool1(
            x, edge_index, edge_attr, batch)
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = self.conv2(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool2(
            x, edge_index, edge_attr, batch)
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = x1 + x2
        x = self.out(x)
        return x
