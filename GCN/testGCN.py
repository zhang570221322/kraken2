from torch.functional import Tensor
from model import GCN
import torch
from torch import nn
# from handle_data import train_data, test_data, load_array


class Arg:
    num_epochs = 100
    lr = learning_rate = 5
    weight_decay = 0
    batch_size = 64


x, adj, y = torch.load(
    "/data/home/wlzhang/classfication/software/kraken2/kraken2_weight/kraken2/GCN/test-data")
train_iter = [[x, adj, y]]
arg = Arg()
net = GCN(adj, 2, 220)
loss = nn.CrossEntropyLoss()

# 保存结果
train_ls = []
real_ls = y.unsqueeze(-1).max(axis=1).indices


def get_prd_label(net, features, labels):
    y_hat = net(features)
    label_hat = y_hat.squeeze().max(axis=1)
    return [label_hat.indices, labels]


# 这里使用的是Adam优化算法
optimizer = torch.optim.Adam(net.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
for epoch in range(arg.num_epochs):
    for X, adj, y in train_iter:
        optimizer.zero_grad()
        l = loss(net((X, adj)), y)
        l.backward()
        optimizer.step()
    train_ls.append(net((X, adj)).squeeze().max(axis=1).indices)
