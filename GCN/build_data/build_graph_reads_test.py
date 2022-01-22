from build_graph_reads import get_feature
import matplotlib.pyplot as plt
import networkx as nx
import torch
import pdb
from torch_geometric.utils import add_self_loops
from torch_geometric.data import data as D
file = "./test.fastq"
temp = get_feature(file)
x, adj, y = next(temp)
u, v, weight = adj
pdb.set_trace()
edge_index = torch.tensor([u, v], dtype=torch.int)   # 2 x E
edge_index = add_self_loops(edge_index, num_nodes=1024)
x = torch.tensor(x, dtype=torch.float).unsqueeze(-1)   # N x emb(in)
edge_attr = torch.tensor(weight, dtype=torch.float)   # E x edge_dim
train_mask = torch.tensor([True], dtype=torch.bool)
val_mask = train_mask
test_mask = train_mask
g = D.Data()
g.x, g.edge_index, g.edge_attr = x, edge_index, edge_attr
fig, ax = plt.subplots()
nx.draw(g.to_networkx(), ax=ax)
ax.set_title('Class: 0')
plt.show()
