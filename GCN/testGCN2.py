
from modelGCN import Net
from common import Arg, my_plot, DataLoaderTqdm
from ReadsDataSet import ReadsDataSet
from sklearn.metrics import roc_auc_score
import torch
from torch import nn
import numpy as np
import pdb

arg = Arg(20, 0.01, 0.75, 1000)
# data
print("start to load data...")
dataset = ReadsDataSet(root='./')
dataset = dataset.shuffle()
len_dataset = len(dataset)
train_dataset = dataset[:len_dataset//10*8]
val_dataset = dataset[len_dataset//10*8:len_dataset//10*9]
test_dataset = dataset[len_dataset//10*9:]
print("load data done!")

train_loader = DataLoaderTqdm(train_dataset, batch_size=arg.batch_size)
val_loader = DataLoaderTqdm(val_dataset, batch_size=arg.batch_size)
test_loader = DataLoaderTqdm(test_dataset, batch_size=arg.batch_size)

input_dim = dataset.num_features
output_dim = dataset.num_classes
# model
# device = torch.device('cuda:1')
device = torch.device('cuda:1')
model = Net(input_dim, output_dim).to(device)
optimizer = torch.optim.Adam(model.parameters(),
                             lr=arg.learning_rate,
                             weight_decay=arg.weight_decay)
crit = torch.nn.CrossEntropyLoss()


def train():
    model.train()
    loss_all = 0
    for data in train_loader:
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        label = data.y.to(device)
        loss = crit(output, label)
        loss.backward()
        loss_all += loss.item()
        optimizer.step()
    return loss_all / len(train_dataset)


def evaluate(loader):
    model.eval()
    predictions = []
    labels = []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            pred = model(data).detach().cpu().numpy()

            label = data.y.detach().cpu().numpy()
            predictions.append(pred)
            labels.append(label)

    predictions = np.hstack(predictions)
    labels = np.hstack(labels)

    return roc_auc_score(labels, predictions)


loss_record = []
train_acc_record = []
val_acc_record = []
test_acc_record = []
for epoch in range(arg.num_epochs):
    loss = train()
    train_acc = evaluate(train_loader)
    val_acc = evaluate(val_loader)
    test_acc = evaluate(test_loader)
    loss_record.append(loss)
    train_acc_record.append(train_acc)
    val_acc_record.append(val_acc)
    test_acc_record.append(loss_record)
    print('Epoch: {:03d}, Loss: {:.5f}, Train Auc: {:.5f}, Val Auc: {:.5f}, Test Auc: {:.5f}'.
          format(epoch, loss, train_acc, val_acc, test_acc))

my_plot(loss_record, train_acc_record,
        val_acc_record, test_acc_record, ["Loss", "Train Auc", "Val Auc", "Test Auc"], arg=arg)
