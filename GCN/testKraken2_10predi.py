#%%
from modelGCN import *
from common import Arg , DataLoaderTqdm 
from ReadsDataSet_kraken_all import ReadsDataSet
import torch,random
import numpy as np
import pdb
from collections import Counter
seed=12345
# device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
device = torch.device('cuda:0')
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)

arg = Arg(100, 0.000942607974725336,0.0009, 1000)
# data
print("start to load data...")
dataset = ReadsDataSet(root='./')
input_dim = dataset.num_features
len_dataset = len(dataset)
dataset_loader = DataLoaderTqdm(dataset, batch_size=arg.batch_size)
print("load data done!")
output_dim = 1
model = Net(input_dim,output_dim).to(device)
model.load_state_dict(torch.load("./args/Kraken2"))




def evaluate(loader):
    model.eval()
    predictions = []
    labels = []
    with torch.no_grad():
        for data in loader:
            data = data.to(device)
            pred = model(data)
            pred = pred.squeeze().round().detach().cpu().numpy()
            label = data.y.detach().cpu().numpy()
            predictions.append(pred)
            labels.append(label)
    predictions = np.concatenate(predictions, axis=0)
    labels = np.concatenate(labels, axis=0)
    acc = (labels==predictions).sum()/len(labels)
    return acc,predictions
print("load model done!")

#%%
acc,predictions = evaluate(dataset_loader)
print("predictions done!")
print(f"acc is {acc}")
#%%

def get_gcnpredi(predictions,dataset):
    
    predi_all =[]
    for g,predi_index in zip(dataset,predictions):
        predi_all.append(g.y_list[int(predi_index)].item())
    index_all = list(dataset.data.seq_index.detach().cpu().numpy())
    
    start =0
    end =0
    predi_file_all=[]
    index_file_all=[]
    for i in range(10):
        
        total = (dataset.data.file_num==i).sum()
        end = total+end
        predi = predi_all[start:end]
        predi_file_all.append(predi)
        index = index_all[start:end]
        index_file_all.append(index)
        start =end
    return index_file_all,predi_file_all
index_file_all,predi_file_all = get_gcnpredi(predictions,dataset)
print("get_gcnpredi to taxid done!")

#%%
RANKS = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
DICT_RANK_TO_INDEX = dict(zip(RANKS, list(range(len(RANKS)))))
RANKS_LOW2HIGH = list(reversed(RANKS))


def load_merged(names_file_path):
    tax_id_to_tax_id = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id_to_tax_id[int(line[0])] = int(line[1])
    return tax_id_to_tax_id


def load_names(tax_id_to_rank, names_file_path):
    tax_id_to_name = {}
    with open(names_file_path) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id = int(line[0])

            if line[3] == "scientific name":
                tax_id_to_name[tax_id] = line[1]
            else:
                continue

            if tax_id_to_rank[tax_id] == "species" or tax_id_to_rank[tax_id] == "strain":
                names = tax_id_to_name[tax_id].split(" ")
                if len(names) > 2:
                    if tax_id_to_rank[tax_id] == "strain":
                        tax_id_to_name[tax_id] = "{} {} strain".format(names[0], names[1])
                    else:
                        tax_id_to_name[tax_id] = "{} {}".format(names[0], names[1])
    tax_id_to_name[""] = ""
    tax_id_to_name[0] = "unassigned"
    return tax_id_to_name


def check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id):
    if tax_id_to_parent[tax_id] in tax_id_to_parent:
        if tax_id_to_rank[tax_id_to_parent[tax_id]] == "species":
            return True
        elif tax_id_to_parent[tax_id] != 1 and tax_id_to_rank[tax_id_to_parent[tax_id]] not in RANKS:
            return check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id_to_parent[tax_id])
    return False


def load_tax_info(ncbi_nodes_file):
    tax_id_to_parent = {}
    tax_id_to_rank = {}
    with open(ncbi_nodes_file) as read_handler:
        for line in read_handler:
            if len(line.strip()) == 0:
                continue
            line = line.split('|')
            line = list(map(str.strip, line))
            tax_id = int(line[0])
            tax_id_to_parent[tax_id] = int(line[1])
            tax_id_to_rank[tax_id] = line[2]

    for tax_id, rank in tax_id_to_rank.items():
        if tax_id_to_rank[tax_id] == "no rank" and check_parent_is_species(tax_id_to_parent, tax_id_to_rank, tax_id):
            tax_id_to_rank[tax_id] = "strain"
    tax_id_to_rank[0]="unassigned"
    return tax_id_to_parent, tax_id_to_rank


def get_id_path(tax_id, tax_id_to_parent, tax_id_to_rank, tax_id_to_tax_id):
    if tax_id ==0:
        return [0]
    if tax_id not in tax_id_to_rank:
        if tax_id_to_tax_id and tax_id in tax_id_to_tax_id:
            tax_id = tax_id_to_tax_id[tax_id]
        else:
            print("Invalid NCBI taxonomic ID: {}".format(tax_id))
            return []

    while tax_id_to_rank[tax_id] not in RANKS:
        tax_id = tax_id_to_parent[tax_id]
        if tax_id == 1:
            return []

    index = DICT_RANK_TO_INDEX[tax_id_to_rank[tax_id]]
    path = [''] * (index + 1)
    path[index] = tax_id

    id = tax_id
    while id in tax_id_to_parent:
        id = tax_id_to_parent[id]
        if id == 1:
            break
        if tax_id_to_rank[id] not in RANKS:
            continue
        index = DICT_RANK_TO_INDEX[tax_id_to_rank[id]]
        path[index] = id
        if tax_id_to_rank[id] == "superkingdom":
            break
    return path
taxonomy_db="/home/yegpu/zwl/taxonomy/"
tax_id_to_tax_id = load_merged(f"{taxonomy_db}/merged.dmp")
tax_id_to_parent,tax_id_to_rank = load_tax_info(f"{taxonomy_db}/nodes.dmp")
tax_id_to_name = load_names(tax_id_to_rank,f"{taxonomy_db}/names.dmp")
taxonomy=[tax_id_to_parent,tax_id_to_rank,tax_id_to_tax_id]
print("load taxonomy done!")
#%%

def restor_data(lines_pool,seq_index_list,predi_list):
    """
    流式输出: (weight,adjacent_matrix,Y)
    """
    gcn_predi_index =0
    for seq_index,line in enumerate(lines_pool):
        line = line.strip()
        line_split = line.split('\t')
        seq_id = line_split[1]
        real_taxid = seq_id.split("|")[-1].split("_")[0]
        predi_taxid = line_split[2]
        predi = 0
        # 排除一个也没命中的
        # try:
        if  gcn_predi_index !=len(seq_index_list) and seq_index == seq_index_list[gcn_predi_index]:
            predi = predi_list[gcn_predi_index]
            predi = real_taxid
            gcn_predi_index+=1
        else:
            predi = real_taxid
        yield seq_id,predi
        # except Exception as e:
        #     print(e)
        #     print(seq_index,gcn_predi_index)
            

from collections import Counter
for type in ["origin","weight"]:
    print(f"handle {type}")
    with open(f"GCN.{type}.OPAL","w") as f:
        for i in range(10):
            print(f"handle {i}")
            version = "1.0.0"
            sample_id=str(i)
            file = f"/home/yegpu/zwl/kraken2/GCN/raw/kraken_output_data_graph/{i}.fastq.{type}.output"
            stac=set()
            for i,j in restor_data(open(file,"r"),index_file_all[i],predi_file_all[i]):
                stac.add(j)
            print(len(stac))
            # f.write(f"@SampleID:{sample_id}\n")
            # f.write(f"@Version:{version}\n")
            # f.write(f"@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n") 
            # f.write(f"\n")
            # f.write(f"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
            # total=0
            # total_taxid=[]
            # for seq_id,predi in  restor_data(open(file,"r"),index_file_all[i],predi_file_all[i]):
            #     taxid=int(predi)
            #     total+=1
            #     taxid_path=get_id_path(taxid,*taxonomy)
            #     total_taxid+=taxid_path
            # temp = dict(Counter(total_taxid))
            # temp.pop('')
            # for taxid,times in temp.items():
            #     rank = tax_id_to_rank[taxid]
            #     taxpath = get_id_path(taxid,*taxonomy)
            #     taxpath_text =  "|".join([str(i) for i in taxpath])
            #     taxpathsn_text = "|".join([str(tax_id_to_name[i]) for i in taxpath])
            #     percentage = str(round(times/total,10))
            #     text = f"{taxid}\t{rank}\t{taxpath_text}\t{taxpathsn_text}\t{percentage}\n"
            #     f.write(text)
# %%
