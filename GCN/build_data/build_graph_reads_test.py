from build_graph_reads import get_feature,DEFAULT_K
import matplotlib.pyplot as plt
import networkx as nx
import torch
import pdb
from torch_geometric.utils import add_self_loops,to_networkx
from torch_geometric.data import data as D
file = "./test.fastq"
temp = get_feature(file)
x, adj, y = next(temp)
u, v, weight = adj
edge_index = torch.tensor([u, v], dtype=torch.int)    
edge_attr = torch.tensor(weight, dtype=torch.float)    
edge_index,edge_attr = add_self_loops(edge_index,edge_attr,fill_value=1, num_nodes=4**DEFAULT_K)
x = torch.tensor(x, dtype=torch.float).unsqueeze(-1)   
train_mask = torch.tensor([True], dtype=torch.bool)
val_mask = train_mask
test_mask = train_mask
g = D.Data()
g.x, g.edge_index, g.edge_attr = x, edge_index, edge_attr
fig, ax = plt.subplots()
nx.draw(to_networkx(g), ax=ax)
ax.set_title(f'Class: {y}')
plt.show()
