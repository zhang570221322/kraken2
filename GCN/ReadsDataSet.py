# %%
from numpy import dtype
from common import date_print
import torch
from torch_geometric.utils import add_self_loops, to_networkx
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
from build_data.build_graph_reads import get_feature, DEFAULT_K, get_read_length
# import matplotlib.pyplot as plt
# import networkx as nx
import pdb


class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        return ["test.fastq"]

    @property
    def processed_file_names(self):
        return ['Reads.dataset']

    def download(self):
        pass

    def process(self):
        g_list = []
        data_list = []
        y_list = []
        for fastq_dir in self.raw_paths:
            reads_length = get_read_length(fastq_dir)
            date_print(f"Process files {fastq_dir}")
            for x, adj, y in tqdm(get_feature(fastq_dir), total=reads_length, unit="read"):
                u, v, weight = adj
                edge_index = torch.tensor([u, v], dtype=torch.int)
                edge_attr = torch.tensor(weight, dtype=torch.float)
                edge_index, edge_attr = add_self_loops(
                    edge_index, edge_attr, fill_value=1, num_nodes=4**DEFAULT_K)
                x = torch.tensor(x, dtype=torch.float).unsqueeze(-1)
                data_list.append((x, edge_index, edge_attr, y))
                y_list.append(y)
        y_list = torch.unique(torch.Tensor(y_list))
        y_label = {int(y_list[i]): i for i in range(
            torch.unique(y_list).numel())}
        date_print(f"Transform data to graph...")
        for x, edge_index, edge_attr, y in tqdm(data_list, unit="graph"):
            y = torch.Tensor([y_label[y]]).long()
            g = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y)
            # fig, ax = plt.subplots()
            # nx.draw(to_networkx(g), ax=ax)
            # ax.set_title(f'Class: {y}')
            # plt.show()
            g_list.append(g)
        data, slices = self.collate(g_list)
        torch.save((data, slices), self.processed_paths[0])


# %%
if __name__ == "__main__":
    dataset = ReadsDataSet("./")

# %%
