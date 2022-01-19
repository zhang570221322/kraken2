from model import GCN, Arg
import torch
from torch import nn
import random
from handle_data import getdata, load_array, my_plot
import pdb
random.seed(32)


arg = Arg(20, 0.01, 0.75, 1000)
# device = torch.device('cuda:1')
device = torch.device('cuda:1')
# data
print("start to load data...")
train_data, test_data = getdata()
# x, adj, y = train_data
# train_count = 10000
# test_count = train_count+1000
# train_iter = load_array(
#     (x[:train_count], adj[:train_count], y[:train_count]), arg.batch_size)
# test_data = (x[train_count:test_count],
#              adj[train_count:test_count], y[train_count:test_count])

train_iter = load_array(train_data, arg.batch_size)
print("load data done!")

# model
net = GCN(18, 18).to(device)
loss = nn.MSELoss()
# 测试数据
test_x, test_adj, test_y = test_data
test_x = test_x.to(device)
test_adj = test_adj.float().to(device)
test_y = test_y.squeeze().to(device)
# Adam
optimizer = torch.optim.Adam(net.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
loss_record = []
acc_record = []
for epoch in range(arg.num_epochs):
    net.train()
    for X, adj, y in train_iter:
        X = X.to(device)
        adj = adj.float().to(device)
        y = y.to(device)
        y=y.squeeze()
        optimizer.zero_grad()
        l = loss(net((X, adj)), y) 
        l.backward()
        optimizer.step()
    net.eval()
    test_rl = net((test_x, test_adj))
    print_loss = float(loss(test_rl, test_y).cpu().detach().numpy())
    real_label = test_rl.max(axis=1).indices.cpu().detach().numpy()
    test_label = test_y.max(axis=1).indices.cpu().detach().numpy()
    acc = (test_label == real_label).sum() / len(test_label)
    loss_record.append(print_loss)
    acc_record.append(acc)
my_plot(list(range(len(loss_record))), loss_record, acc_record, arg=arg)
