# %%
from common import date_print, standardization
import torch
from torch_geometric.utils import add_self_loops, to_networkx
from torch_geometric.data import Data, InMemoryDataset
from tqdm import tqdm
from build_data.build_graph_reads import get_feature, get_feature3, Y_Handle
from build_data.graph_model import ReadGenerator, DEFAULT_K
from collections import Counter
import pickle
import matplotlib.pyplot as plt
import networkx as nx
import pdb
file_name_prefix = "10class"
#%%

class ReadsDataSet(InMemoryDataset):
    def __init__(self, root,  transform=None, pre_transform=None):
        super(ReadsDataSet, self).__init__(
            root, transform, pre_transform)
        self.data, self.slices = torch.load(self.processed_paths[0])
        with open(self.processed_paths[0]+"_level_tax", 'rb') as f:
            temp = pickle.load(f)
            self.level_tax = temp[0]
            self.level_dim = temp[1]

    @property
    def raw_file_names(self):
        return [f"{file_name_prefix}.fasta"]

    @property
    def processed_file_names(self):
        return [f'{file_name_prefix}_k{DEFAULT_K}_fasta_dbg.dataset']

    def download(self):
        pass

    def process(self):
        g_list = []
        print(self.__class__)
        for fasta_dir in self.raw_paths:
            reads_length = ReadGenerator(fasta_dir, "fasta2fastq").get_fasta2read_length()
            y_handle = Y_Handle(fasta_dir,"fasta")
            y_handle.muti_label_mode(mode=2)
            print(y_handle.ys_dic)
            print(y_handle.tree)
            date_print(
                f"Process files {fasta_dir}, Transform data to graph...")
            for x, adj, y,_ in tqdm(get_feature3(fasta_dir), total=reads_length, unit="read"):
                # y = get_genus(y)
                # if y in [1870884,1485]:
                #     continue
                y_index,ncbi_lineage = y_handle.get_linear_one_hot(y)
                tax=torch.tensor(ncbi_lineage)
                u, v, weight = adj
                edge_index, edge_attr = torch.tensor(
                    [u, v]),  torch.tensor(weight)
                # edge_index = edge_index.short()
                # edge_attr = edge_attr.short()
                x = torch.tensor(x).flatten(1)
                y = torch.tensor(y_index)
                g = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y,tax1=tax[-1],tax2=tax[-2],tax3=tax[-3])
                g_list.append(g)
                # pdb.set_trace()
                # # #vision
                # plt.figure()
                # G = to_networkx(g)
                # pos = nx.spring_layout(G)
                # nx.draw(G,pos=pos,node_size=30)
                # nx.draw_networkx_labels(G,pos = pos, labels = {})
                # plt.show()
        data, slices = self.collate(g_list)
        print(Counter(data.tax1.numpy()))
        torch.save((data, slices), self.processed_paths[0])
        with open(self.processed_paths[0]+"_level_tax", 'wb') as f:
            pickle.dump([y_handle.level_tax,y_handle.level_dim], f)


# %%
if __name__ == "__main__":
    from common import get_tax_list 
    fastq_dir = "/home/yegpu/zwl/kraken2/GCN/raw/10class.fasta"
    y_handle = Y_Handle(fastq_dir,"fasta")
    y_handle.muti_label_mode(mode=2)
    dataset = ReadsDataSet("./")
    taxs=dataset.level_tax
    error=[]
    for data in dataset:
        try:
            temp=get_tax_list(taxs,list(data.y.numpy()),mode=2)
            
        except KeyError as err:
                pdb.set_trace()
                raise err    
        data_tax = [data.tax3,data.tax2,data.tax1]
        for i in range(1,4):
            real= int(temp[-i])
            if real ==0:
                    continue
            predi=data_tax.pop()
            if real !=  predi.item():
                error.append(data)
    print(len(error))
# %%
