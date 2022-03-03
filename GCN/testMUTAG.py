
from torch_geometric.loader import DataLoader
from torch_geometric.datasets import TUDataset
from modelGCN import *
import random
import numpy as np
from common import my_plot
from collections import Counter
seed = 12345
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)

dataset = TUDataset(root='data/TUDataset', name='MUTAG')
dataset = dataset.shuffle()

train_dataset = dataset[:150]
test_dataset = dataset[150:]
train_loader = DataLoader(train_dataset, batch_size=64, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)
input_dim = dataset.num_features
output_dim = dataset.num_classes
model = Net2(input_dim, output_dim).to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
criterion = torch.nn.CrossEntropyLoss()


def train():
    model.train()
    loss_all = 0
    for data in train_loader:  # Iterate in batches over the training dataset.
        data = data.to(device)
        out = model(data)  # Perform a single forward pass.
        loss = criterion(out, data.y)  # Compute the loss.
        loss.backward()  # Derive gradients.
        loss_all += loss.item()
        optimizer.step()  # Update parameters based on gradients.
        optimizer.zero_grad()  # Clear gradients.
    return loss_all / len(train_dataset)


def evaluate(loader):
    model.eval()
    predictions = []
    labels = []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            out = model(data)
            pred = out.detach().cpu().numpy()
            label = data.y.detach().cpu().numpy()
            predictions.append(pred)
            labels.append(label)
    predictions = np.concatenate(predictions, axis=0)
    predictions = predictions.argmax(axis=-1)
    print(Counter(predictions), end=",")
    labels = np.concatenate(labels, axis=0)
    print(Counter(labels), end=",")
    acc = (labels == predictions).sum()/len(labels)
    return acc


loss_record = []
train_acc_record = []
test_acc_record = []
for epoch in range(1, 15):
    loss = train()
    test_acc = evaluate(test_loader)
    loss_record.append(loss)
    test_acc_record.append(test_acc)
    message = f'Epoch: {epoch+1}, Loss: {loss:.5f},   Test Auc: {test_acc:.5f}'
    print(message)
my_plot(loss_record, test_acc_record, ["Loss",  "Test Acc"])
