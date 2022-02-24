# %%
from common import date_print , standardization
import torch
from torch_geometric.utils import add_self_loops, to_networkx
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
from build_data.build_graph_reads import get_feature,get_feature2, DEFAULT_K, get_read_length
from ete3 import NCBITaxa
from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
import pdb
file_name_prefix="5class"
ncbi = NCBITaxa()
taxid_levelup_taxid={}
def get_genus(tax_id):
    target = taxid_levelup_taxid.get(tax_id,None)
    if target!=None:
        return target
    loc = -1
    temp = ncbi.get_lineage(tax_id)
    target = temp[loc]
    while ncbi.get_rank([target])[target] != "genus":
        loc = loc-1
        target = temp[loc]
    taxid_levelup_taxid[tax_id]=target
    return target

y_dic = {-1:0}
def one_hot_index(y):
    y_index = y_dic.get(y)
    if y_index == None:
        y_dic[y] = y_dic[-1]
        y_index = y_dic[-1]
        y_dic[-1] =y_dic[-1]+1
    return y_index
class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])
    @property
    def raw_file_names(self):
        return [f"{file_name_prefix}.fastq"]

    @property
    def processed_file_names(self):
        return [f'{file_name_prefix}_k{DEFAULT_K}_all.dataset']

    def download(self):
        pass

    def process(self):
        g_list = []
        print(self.__class__)
        for fastq_dir in self.raw_paths:
            reads_length = get_read_length(fastq_dir)
            date_print(f"Process files {fastq_dir}, Transform data to graph...")
            for x, adj, y,node_labels in tqdm(get_feature2(fastq_dir), total=reads_length, unit="read"):
                # y = get_genus(y)
                # if y in [1870884,1485]:
                #     continue
                y_index = one_hot_index(y)
                u, v, weight = adj
                edge_index, edge_attr =  torch.tensor([u, v]) ,  torch.tensor(weight)
                edge_index=edge_index.long()
                edge_attr = edge_attr.short()
                x = torch.tensor(x, dtype=torch.float)
                # x = standardization(x)
                y = torch.tensor([y_index])
                g= Data(x=x, edge_index=edge_index, edge_attr=None, y=y)
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
        print(Counter(data.y.numpy()))
        print(y_dic)
        torch.save((data, slices), self.processed_paths[0])


# %%
if __name__ == "__main__":
    dataset = ReadsDataSet("./")

# %%