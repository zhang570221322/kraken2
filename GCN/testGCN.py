from model_line import FC, Arg
import torch
from torch import nn
import random
from handle_data import getdata, load_array, my_plot,de_to_dense
import pdb
import numpy as np
random.seed(32)


arg = Arg(30, 0.005, 0.75, 200)
# device = torch.device('cpu')
device = torch.device('cuda:1')
# data
print("start to load data...")
train_data,output_dim, test_data = getdata()
 

train_iter = load_array(train_data, arg.batch_size)
print("load data done!")
# input_dim = train_data[0].shape[2]
# model
net = FC(1024, output_dim).to(device)
# loss = nn.CrossEntropyLoss()
loss = nn.MSELoss()
# 测试数据
test_x,test_adj,test_y=test_data
# test_iter = load_array(test_data, arg.batch_size*2)

def test_net(X,adj,test_y):
    X = X.float().to(device)
    adj = adj.float().to(device)
    test_y =test_y.to(device)
    predi_rl=net((X, adj))
    l = loss(predi_rl, test_y) 
    print_loss = float(l)
    test_y = test_y.cpu().detach().numpy()
    predi_label = predi_rl.max(axis=1).indices.cpu().detach().numpy()
    acc = (predi_label == test_y).sum() / len(test_y)
    return print_loss,acc
# Adam
optimizer = torch.optim.Adam(net.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
loss_record = []
acc_record = []
for epoch in range(arg.num_epochs):
    net.train()
    for X, adj, y in train_iter:
        # adj=de_to_dense(adj)
        X = X.to(device).float()
        adj = adj.to(device).float()
        y = y.to(device)
        res=net((X, adj))
        l = loss(res, y) 
        optimizer.zero_grad()
        l.backward()
        optimizer.step()
    net.eval()
    # pdb.set_trace()
    # for X, adj, test_y in test_iter:
    print_loss,acc =test_net(test_x,test_adj,test_y)
    loss_record.append(print_loss)
    acc_record.append(acc)
    print(print_loss," ",acc)
my_plot(list(range(len(loss_record))), loss_record, acc_record, arg=arg)
