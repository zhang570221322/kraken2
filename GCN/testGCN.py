from model import GCN
import torch
from torch import nn
import random
from handle_data import my_plot
import pdb
random.seed(32)


class Arg:
    num_epochs = 200
    lr = learning_rate = 5
    weight_decay = 0
    batch_size = 256


arg = Arg()
# device = torch.device('cuda:1')
device = torch.device('cpu')
# data
print("start to load data...")
x, adj, y = torch.load("./test-data")
y = y.unsqueeze(-1)
train_data = test_data = (x, adj, y)
train_iter = [[x, adj, y]]
print("load data done!")
# model
net = GCN(2, 1).to(device)
loss = nn.SmoothL1Loss()
# 测试数据
test_x, test_adj, test_y = test_data
test_x = test_x.to(device)
test_adj = test_adj.float().to(device)
test_y = test_y.to(device)
real_ls = test_y.unsqueeze(-1).max(axis=1).indices.T
# Adam
optimizer = torch.optim.Adam(net.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
loss_record = []
acc_record = []
for epoch in range(arg.num_epochs):
    for X, adj, y in train_iter:
        X = X.to(device)
        adj = adj.float().to(device)
        y = y.to(device)
        optimizer.zero_grad()
        l = loss(net((X, adj)), y)
        l.backward()
        optimizer.step()
    print_loss = float(loss(net((test_x, test_adj)),
                       test_y).cpu().detach().numpy())
    real_label = net((test_x, test_adj)).squeeze().max(
        axis=1).indices.cpu().detach().numpy()
    test_label = test_y.squeeze().max(axis=1).indices.cpu().detach().numpy()
    acc = (test_label == real_label).sum() / len(test_label)
    loss_record.append(print_loss)
    acc_record.append(acc)
my_plot(list(range(len(loss_record))), loss_record, acc_record)
