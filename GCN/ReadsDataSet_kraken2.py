# %%
from common import date_print  
import torch
from torch_geometric.utils import add_self_loops, to_networkx
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
import matplotlib.pyplot as plt
import networkx as nx,numpy as np
import pdb



class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])
    @property
    def raw_file_names(self):
        return [f"1.fastq"]

    @property
    def processed_file_names(self):
        return ['kraken2.dataset']

    def download(self):
        pass

    def process(self):
        g_list = []
        preifx="/home/yegpu/zwl/data/kraken2_data_gcn"
        x1 = np.load(f"{preifx}/x1.npy", allow_pickle=True)
        # x2 = np.load(f"{preifx}/x2.npy", allow_pickle=True)
        adjacent_matrixs = np.load( f"{preifx}/adjacent_matrixs.npy", allow_pickle=True)
        Y = np.load(f"{preifx}/Y.npy", allow_pickle=True)
        date_print(f"Process files {preifx}, Transform data to graph...")
        for x, adj, y in tqdm(zip(x1,adjacent_matrixs,Y), total=len(x1), unit="read"):
            adj = np.array(adj)
            u, v = list(adj[:,0]),list(adj[:,1])
            edge_index ,_= add_self_loops(
                    torch.tensor([u, v]), fill_value=1, num_nodes=len(x))
            x = torch.tensor(x, dtype=torch.float).unsqueeze(-1)
            y = torch.tensor(y).argmax()
            g= Data(x=x, edge_index=edge_index, y=y)
            g_list.append(g)
            # pdb.set_trace()
            # vision
            # plt.figure()
            # G = to_networkx(g)
            # pos = nx.spring_layout(G)
            # nx.draw(G,pos=pos,node_size=30)
            # nx.draw_networkx_labels(G,pos = pos, labels = {}) 
            # plt.show()
        data, slices = self.collate(g_list)
        torch.save((data, slices), self.processed_paths[0])


# %%
if __name__ == "__main__":
    dataset = ReadsDataSet("./")

# %%
