# %%
from common import date_print
import torch
from torch_geometric.utils import add_self_loops, to_networkx
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
from build_data.build_graph_kranken2_output import get_feature,get_total_tqdm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pdb
import sys
# file_number = sys.argv[1]
file_number = "0"
print(file_number)

class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        
        # return  [f"{i}.fastq.origin.output" for i in range(10)]
        return  [f"{file_number}.fastq.origin.output"]
    @property
    def processed_file_names(self):
        return  [f"{file_number}.kraken_output_data"]

    def download(self):
        pass

    def process(self):
        g_list = []
        for fastq_dir in self.raw_paths:
            date_print(f"{self.__class__}  Process files {fastq_dir}, Transform data to graph...")
            total = get_total_tqdm(fastq_dir)
            for x, u,v, y ,y_list,seq_index in tqdm(get_feature(fastq_dir), total=total, unit="read"):
                edge_index, _ = add_self_loops(
                    torch.tensor([u, v]), fill_value=1, num_nodes=len(x))
                x = torch.tensor(x, dtype=torch.float).unsqueeze(-1)
                y = torch.tensor(y)
                y_list =  torch.tensor(y_list)
                seq_index = torch.tensor(seq_index)
                g = Data(x=x, edge_index=edge_index, y=y, y_list=y_list,seq_index=seq_index)
                g_list.append(g)
                # pdb.set_trace()
                # # vision
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
