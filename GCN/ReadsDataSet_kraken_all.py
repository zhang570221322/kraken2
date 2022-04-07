# %%
from common import date_print
import torch
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pdb
import sys
class ReadsDataSet_split(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet_split, self).__init__(
            root[0], transform, pre_transform)
        self.data, self.slices = torch.load(root[1])

    @property
    def raw_file_names(self):
        return  ["0.kraken_output_data"]
    @property
    def processed_file_names(self):
        return  [f"0.kraken_output_data"]

    def download(self):
        pass
    def process(self):
        pass
    
class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return  [f"{i}.kraken_output_data" for i in range(10)]
    @property
    def processed_file_names(self):
        return  [f"all.kraken_output_data"]

    def download(self):
        pass

    def process(self):
        g_list = []
        for file_num,fastq_dir in enumerate(self.raw_paths):
            date_print(f"{self.__class__}  Process files {fastq_dir}, Transform data to graph...")
            dataset =  ReadsDataSet_split(["./",fastq_dir])
            
            for g in tqdm(dataset, total=len(dataset), unit="read"):
                g.file_num =torch.tensor(file_num)
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
