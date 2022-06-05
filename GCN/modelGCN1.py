from re import S
import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.nn import  GraphConv, GCNConv, TopKPooling
from torch_geometric.nn import global_mean_pool as gap, global_max_pool as gmp
import pdb

GCN = GCNConv


class Net(nn.Module):
    def __init__(self, input_dim, output_dim) -> None:
        super().__init__()
        hcs = 512
        self.output_dim= output_dim
        self.conv1 = GCN(input_dim, hcs)
        self.pool1 = TopKPooling(hcs, ratio=0.9)
        self.conv2 = GCN(hcs, hcs)
        self.pool2 = TopKPooling(hcs, ratio=0.9)
        self.conv3 = GCN(hcs, hcs)
        self.pool3 = TopKPooling(hcs, ratio=0.9)
        self.conv3 = GCN(hcs, hcs)
        self.pool3 = TopKPooling(hcs, ratio=0.9)
        self.conv4 = GCN(hcs, hcs)
        self.pool4 = TopKPooling(hcs, ratio=0.9)
        self.conv5 = GCN(hcs, hcs)
        self.pool5 = TopKPooling(hcs, ratio=0.9)
        self.conv6 = GCN(hcs, hcs)
        self.pool6 = TopKPooling(hcs, ratio=0.9)
        self.conv7 = GCN(hcs, hcs)
        self.pool7 = TopKPooling(hcs, ratio=0.9)
        self.norm = nn.LayerNorm(input_dim)
        i=0
        self.out0=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out1=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out2=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out3=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out4=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out5=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        i+=1
        self.out6=nn.Sequential(
            nn.ReLU(),
            nn.Linear(hcs*2, hcs*8), nn.ReLU(),
            nn.Linear(
                hcs*8, hcs*32), nn.ReLU(),
            nn.Linear(hcs*32, output_dim[i]) )
        if len(output_dim)==8:
            i+=1
            self.out7=nn.Sequential(
                nn.ReLU(),
                nn.Linear(hcs*2, hcs*8), nn.ReLU(),
                nn.Linear(
                    hcs*8, hcs*32), nn.ReLU(),
                nn.Linear(hcs*32, output_dim[i]) )

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        edge_index = edge_index.long()
        x=x.float()
        if edge_attr != None:
            edge_attr = edge_attr.float()
        edge_attr = None
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
        
        x = self.conv4(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool4(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x4 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv5(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool5(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x5 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv6(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool6(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x6 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        
        x = self.conv7(x, edge_index, edge_attr)
        x = x.relu()
        x, edge_index, edge_attr, batch, _, _ = self.pool7(
            x, edge_index, edge_attr, batch)
        # pdb.set_trace()
        x7 = torch.cat([gmp(x, batch), gap(x, batch)], dim=1)
        x = x1 + x2 + x3 + x4 +x5 + x6+ x7
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
        if len(self.output_dim)==8:
            i+=1
            res.append(self.out7(x))
        return res
    def criterion(self,loss_func,outputs,tar,weight=None):
        losses = 0
        if weight==None:
            weight = torch.ones((len(outputs,)))
        for i  in range(len(outputs)):
            losses+=(loss_func(outputs[i], tar[:,i])*weight[i])
        return losses
