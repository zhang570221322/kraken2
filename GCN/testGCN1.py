# %%
from modelGCN1 import *
from autogl.datasets import utils
from torch_geometric.loader import DataLoader
from common import Arg, my_plot, DataLoaderTqdm, log_time, get_tax_list
from ReadsDataSet1 import ReadsDataSet
import torch
import random
from torch import nn
import numpy as np
import pdb
from collections import Counter

seed = 12345
# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
device = torch.device('cuda:1')
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)

arg = Arg(30, 0.000942607974725336, 0.75, 61)
# data
print("start to load data...")
dataset = ReadsDataSet(root='./')
dataset = utils.graph_random_splits(
    dataset, train_ratio=0.8, val_ratio=0.1, seed=seed)
len_dataset = len(dataset)
train_dataset = dataset.train_split
test_dataset = dataset.test_split
# %%
train_loader = DataLoaderTqdm(train_dataset, batch_size=arg.batch_size)
test_loader = DataLoaderTqdm(test_dataset, batch_size=arg.batch_size)
print("load data done!")
input_dim = dataset.num_features
output_dim = 11
# pdb.set_trace()
model = Net(input_dim, output_dim).to(device)
optimizer = torch.optim.Adam(model.parameters(),
                             lr=arg.learning_rate)
#  weight_decay=arg.weight_decay)
crit = nn.BCELoss()
# @log_time("train")


def train():
    model.train()
    loss_all = 0
    for data in train_loader:
        data = data.to(device)
        output = model(data)
        label = data.y.to(device)
        reshape_len = int(len(label)/output_dim)
        label = label.reshape(reshape_len, output_dim).float()
        # pdb.set_trace()
        loss = crit(output, label)
        loss.backward()
        loss_all += loss.item()
        optimizer.step()
        optimizer.zero_grad()
        # print(loss.item())
    return loss_all / len(train_dataset)

# @log_time("evaluate")


def evaluate(loader):
    model.eval()
    predictions = []
    labels = []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            pred = model(data).round().detach().cpu().numpy()
            label = data.y.detach().cpu().numpy()
            predictions.append(pred)
            labels.append(label)
    predictions = np.concatenate(predictions, axis=0).flatten()
    labels = np.concatenate(labels, axis=0)
    acc = (labels == predictions).sum()/len(labels)
    return acc
# def get_right_tax(data):


loss_record = []
train_acc_record = []
test_acc_record = []
for epoch in range(arg.num_epochs):
    loss = train()
    train_acc = evaluate(train_loader)
    test_acc = evaluate(test_loader)
    loss_record.append(loss)
    train_acc_record.append(train_acc)
    test_acc_record.append(test_acc)
    message = f'Epoch: {epoch+1}, Loss: {loss:.5f}, Train Auc: {train_acc:.5f} , Test Auc: {test_acc:.5f}'
    print(message)
my_plot(loss_record, train_acc_record, test_acc_record,
        ["Loss", "Train Acc", "Test Acc"], arg=arg)
torch.save(model.state_dict(), "./args/Bruijn_k28")
