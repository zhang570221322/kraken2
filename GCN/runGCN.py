from model import GCN
import torch
from torch import nn
import random
from handle_data import getdata, load_array
random.seed(32)


class Arg:
    num_epochs = 30
    lr = learning_rate = 5
    weight_decay = 0
    batch_size = 256


arg = Arg()
device = torch.device('cuda:1')

# data
print("start to load data...")
train_data, test_data = getdata()
train_iter = load_array(train_data, arg.batch_size)
print("load data done!")
# model
net = GCN(2, 1).to(device)
loss = nn.SmoothL1Loss()
# 测试数据
test_x, test_adj, test_y = test_data
test_x = test_x.to(device)
test_adj = test_adj.to(device)
test_y = test_y.to(device)
real_ls = test_y.unsqueeze(-1).max(axis=1).indices.T
# Adam
optimizer = torch.optim.Adam(net.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
for epoch in range(arg.num_epochs):
    for X, adj, y in train_iter:
        X = X.to(device)
        adj = adj.float().to(device)
        y = y.to(device)
        optimizer.zero_grad()
        l = loss(net((X, adj)), y)
        l.backward()
        optimizer.step()
    print_loss = loss(net((test_x, test_adj)), test_y.unsqueeze(-1))
    print(f"epoch-{epoch}:loss is {print_loss}")
    train_ls = net((test_x, test_adj)).squeeze().max(axis=1).indices
    print(train_ls)
