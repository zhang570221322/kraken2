from model import GCN, Arg
import torch
from torch import nn
import random
from handle_data import getdata, load_array, my_plot
import pdb
random.seed(32)


args = [Arg(60, 0.01, 0.88, 14240),
        # Arg(200, 0.05, 0.78, 512),
        # Arg(400, 0.01, 0.98, 256),
        # Arg(400, 0.01, 0.75, 256),
        # Arg(800, 0.005, 0.78, 512),
        # Arg(800, 0.005,  0.98, 256),
        # Arg(1000, 0.001, 0.75, 512),
        # Arg(1000, 0.001, 0.88, 256),
        ]


device = torch.device('cuda:1')

# data
print("start to load data...")
train_data, test_data = getdata()

print("load data done!")
# model
net = GCN(1, 1).to(device)
loss = nn.MSELoss()
# 测试数据
test_x, test_adj, test_y = test_data
test_x = test_x.to(device)
test_adj = test_adj.float().to(device)
test_y = test_y.to(device)
# Adam
for arg in args:
    train_iter = load_array(train_data, arg.batch_size)
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
        # pdb.set_trace()
        print(epoch)
        print(print_loss)
        print(acc)
        loss_record.append(print_loss)
        acc_record.append(acc)
    my_plot(list(range(len(loss_record))), loss_record, acc_record, arg=arg)
