import os
import pdb
from common import log_time
from build_data.graph_model import feature_space_init, ReadGenerator, DEFAULT_K
import copy
import sys


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
    read_kmers = []
    u, v, weight = [], [], []
    u_v_weight_dic = {}
    for read in read_generator.read_Generator():
        tax_id = 0
        feature_space = copy.deepcopy(feature_space_init)
        tax_id = read.id.split("|kraken:taxid|")[-1].split("_")[0]
        for kmer, _ in read_generator.Kmer_index_Generator(read):
            feature_space[kmer][1] += 1
            read_kmers.append(kmer)
        for kmer, kmer_read_index in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_index = feature_space[kmer][0]
            for i in range(1, len(read.seq)):
                if start-i >= 0:
                    fea_tar_index = feature_space[read_kmers[start-i]][0]
                    key = f"{feature_index},{fea_tar_index}"
                    value = u_v_weight_dic.get(key)
                    if not value or -i-1 > value:
                        u.append(feature_index)
                        v.append(fea_tar_index)
                        weight.append(-i-1)
                        u_v_weight_dic[key] = -i-1
                if end + i <= len(read.seq):
                    fea_tar_index = feature_space[read_kmers[start+i]][0]
                    key = f"{feature_index},{fea_tar_index}"
                    value = u_v_weight_dic.get(key)
                    if not value or i+1 > value:
                        u.append(feature_index)
                        v.append(fea_tar_index)
                        weight.append(i+1)
                        u_v_weight_dic[key] = i+1
        if tax_id == 0 or (not tax_id.isdigit()):
            continue
        x = [temp[1] for temp in list(feature_space.values())]
        yield x, (u, v, weight), int(tax_id)
        read_kmers.clear()
        u.clear()
        v.clear()
        weight.clear()
        u_v_weight_dic.clear()


def get_read_length(file):
    temp = os.popen(f"wc -l {file}").readline().split(" ")[0]
    if temp.isdigit():
        return int(temp)//4
    else:
        raise Exception(f"wc -l {file} return error")
