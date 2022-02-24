import os
import pdb,sys
from common import log_time,AUTO_INCREMENT
from build_data.graph_model import feature_space_init, ReadGenerator, DEFAULT_K,feature_str
import copy
from common import standardization
kmer_encode_dic={}

def kmer_encode(kmer:str):
    temp = kmer_encode_dic.get(kmer,None)
    if temp==None:
        feature = feature_str(kmer)
        kmer_encode_dic[kmer]=feature
        return feature
    else:
        return temp
    
def data_iter(batch, data):
    """
    # 均分data
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]


# def set_value(adj_matrix, feature_index, read_kmers, end, tar_value, feature_space):
#     fea_target_feature = feature_space[read_kmers[end]][0]
#     value = abs(adj_matrix[feature_index][fea_target_feature])
#     if value == 0:
#         adj_matrix[feature_index][fea_target_feature] = tar_value
#     elif abs(value) > abs(tar_value):
#         adj_matrix[feature_index][fea_target_feature] = tar_value
 

def get_feature(file_name):
    read_generator = ReadGenerator(file_name, "reads")
    feature_space_init.init()
    read_kmers = []
    u, v, weight = [], [], []
    u_v_weight_dic = {}
    for read in read_generator.read_Generator():
        tax_id = 0
        feature_space = copy.deepcopy(feature_space_init)
        tax_id = read.id.split("|kraken:taxid|")[-1].split("_")[0]
        for kmer, kmer_read_index,_ in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_space[kmer][1] += 1
            feature_space[kmer][2] = start
            read_kmers.append(kmer)
        for kmer, kmer_read_index,_ in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_index = feature_space[kmer][0]
            # for i in range(1,len(read.seq)):
            for i in range(1,2):
                if start-i >= 0:
                    fea_tar_index = feature_space[read_kmers[start-i]][0]
                    key = f"{feature_index},{fea_tar_index}"
                    value = u_v_weight_dic.get(key)
                    tar_value= i
                    if not value or tar_value > value:
                        u.append(feature_index)
                        v.append(fea_tar_index)
                        weight.append(tar_value)
                        u_v_weight_dic[key] = tar_value
                if end + i <= len(read.seq):
                    fea_tar_index = feature_space[read_kmers[start+i]][0]
                    key = f"{feature_index},{fea_tar_index}"
                    value = u_v_weight_dic.get(key)
                    tar_value = i
                    if not value or tar_value > value:
                        u.append(feature_index)
                        v.append(fea_tar_index)
                        weight.append(tar_value)
                        u_v_weight_dic[key] = tar_value
        if tax_id == 0 or (not tax_id.isdigit()):
            continue
        x = [temp[1:] for temp in list(feature_space.values())]
        yield x, (u, v, weight), int(tax_id),read_kmers
        read_kmers.clear()
        u.clear()
        v.clear()
        weight.clear()
        u_v_weight_dic.clear()



def get_feature2(file_name):
    read_generator = ReadGenerator(file_name, "reads")
    x,u, v, weight =[], [], [], []
    kmer_loc = {}
    # loc_kmer={}
    for read in read_generator.read_Generator():
        tax_id = 0
        pre_kmer=""
        auto_increse=AUTO_INCREMENT()
        tax_id = read.id.split("|kraken:taxid|")[-1].split("_")[0]
        for cur_kmer, kmer_read_index,cur_kmer_feature in read_generator.Kmer_index_Generator(read):
            cur_loc = kmer_loc.get(cur_kmer)
            pre_loc = kmer_loc.get(pre_kmer)
            if cur_loc==None:
                index = next(auto_increse)
                kmer_loc[cur_kmer]=index
                cur_loc = index
                # loc_kmer[start]=cur_kmer
                x.append(cur_kmer_feature)
            # kmer muti 
            # pdb.set_trace()
            if pre_loc != None:
                u.append(pre_loc)
                v.append(cur_loc)
            pre_kmer=cur_kmer
        if tax_id == 0 or (not tax_id.isdigit()):
            continue
        yield x, (u, v, weight), int(tax_id), {v:k for k,v in kmer_loc.items()}
        u.clear()
        v.clear()
        weight.clear()
        x.clear()
        kmer_loc.clear()


def get_read_length(file):
    temp = os.popen(f"wc -l {file}").readline().split(" ")[0]
    if temp.isdigit():
        return int(temp)//4
    else:
        raise Exception(f"wc -l {file} return error")

# %%
