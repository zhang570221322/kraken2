from model import GCN, Arg
import torch
from torch import nn
import random
from handle_data import getdata, load_array, my_plot
import pdb
random.seed(32)


arg = Arg(200, 0.01, 0.88, 256)
# device = torch.device('cuda:1')
device = torch.device('cuda:0')
# data
print("start to load data...")
train_data,_  = getdata()
x, adj, y = train_data
train_iter = [[x[:800], adj[:800], y[:800]]]
test_data = (x[800:1000], adj[800:1000], y[800:1000])
pdb.set_trace()
print("load data done!")
# model
net = GCN(2, 1).to(device)
loss = nn.MSELoss()
# 测试数据
test_x, test_adj, test_y = test_data
test_x = test_x.to(device)
test_adj = test_adj.float().to(device)
test_y = test_y.to(device)
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
        y = y.squeeze()
        l = loss(net((X, adj)), y)
        l.backward()
        optimizer.step()
    # pdb.set_trace()
    print_loss = float(loss(net((test_x, test_adj)),
                       test_y).cpu().detach().numpy())
    real_label = net((test_x, test_adj)).squeeze().max(
        axis=1).indices.cpu().detach().numpy()
    test_label = test_y.squeeze().max(axis=1).indices.cpu().detach().numpy()
    acc = (test_label == real_label).sum() / len(test_label)
    loss_record.append(print_loss)
    acc_record.append(acc)
my_plot(list(range(len(loss_record))), loss_record, acc_record, arg=arg)
