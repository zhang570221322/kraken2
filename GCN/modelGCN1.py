from re import S
import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import GINConv, GraphConv, GCNConv, TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
import pdb

GCN = GraphConv


class Net(nn.Module):
    def __init__(self, input_dim, output_dim) -> None:
        super().__init__()
        hcs = 30
        self.output_dim= output_dim
        self.conv1 = GCN(input_dim, hcs)
        self.pool1 = TopKPooling(hcs, ratio=0.9)
        self.conv2 = GCN(hcs, hcs)
        self.pool2 = TopKPooling(hcs, ratio=0.9)
        self.conv3 = GCN(hcs, hcs)
        self.pool3 = TopKPooling(hcs, ratio=0.9)
        self.norm = nn.LayerNorm(input_dim)
        i=0
        self.out0=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out1=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out2=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out3=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out4=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out5=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out6=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out7=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out8=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out9=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out10=nn.Sequential(
            nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(), nn.Dropout(p=0.15074803632216038),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        edge_index = edge_index.long()
        if edge_attr != None:
            edge_attr = edge_attr.float()
        # pdb.set_trace()
        x = self.norm(x)
        x = self.conv1(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool1(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x1 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = self.conv2(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool2(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x2 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = self.conv3(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool3(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x3 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = x1 + x2 + x3
        res=[]
        i=0
        res.append(self.out0(x))
        i+=1
        res.append(self.out1(x))
        i+=1
        res.append(self.out2(x))
        i+=1
        res.append(self.out3(x))
        i+=1
        res.append(self.out4(x))
        i+=1
        res.append(self.out5(x))
        i+=1
        res.append(self.out6(x))
        i+=1
        res.append(self.out7(x))
        i+=1
        res.append(self.out8(x))
        i+=1
        res.append(self.out9(x))
        i+=1
        res.append(self.out10(x))
        return res
    def criterion(self,loss_func,outputs,tar,weight=None):
        losses = 0
        if weight==None:
            weight = torch.ones((len(outputs,)))
        for i  in range(len(outputs)):
            losses+=(loss_func(outputs[i], tar[:,i])*weight[i])
        return losses
