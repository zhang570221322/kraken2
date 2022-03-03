
import pdb
from tkinter.tix import Tree
from common import log_time, AUTO_INCREMENT
from build_data.graph_model import feature_space_init, ReadGenerator,  feature_str
import copy
from ete3 import NCBITaxa
import numpy as np
ncbi = NCBITaxa()

kmer_encode_dic = {}


def kmer_encode(kmer: str):
    temp = kmer_encode_dic.get(kmer, None)
    if temp == None:
        feature = feature_str(kmer)
        kmer_encode_dic[kmer] = feature
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
        for kmer, kmer_read_index, _ in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_space[kmer][1] += 1
            feature_space[kmer][2] = start
            read_kmers.append(kmer)
        for kmer, kmer_read_index, _ in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_index = feature_space[kmer][0]
            # for i in range(1,len(read.seq)):
            for i in range(1, 2):
                if start-i >= 0:
                    fea_tar_index = feature_space[read_kmers[start-i]][0]
                    key = f"{feature_index},{fea_tar_index}"
                    value = u_v_weight_dic.get(key)
                    tar_value = i
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
        yield x, (u, v, weight), int(tax_id), read_kmers
        read_kmers.clear()
        u.clear()
        v.clear()
        weight.clear()
        u_v_weight_dic.clear()


def get_feature2(file_name):
    read_generator = ReadGenerator(file_name, "reads")
    x, u, v, weight = [], [], [], []
    kmer_loc = {}
    # loc_kmer={}
    for read in read_generator.read_Generator():
        tax_id = 0
        pre_kmer = ""
        auto_increse = AUTO_INCREMENT()
        tax_id = read.id.split("|kraken:taxid|")[-1].split("_")[0]
        for cur_kmer, kmer_read_index, cur_kmer_feature in read_generator.Kmer_index_Generator(read):
            cur_loc = kmer_loc.get(cur_kmer)
            pre_loc = kmer_loc.get(pre_kmer)
            if cur_loc == None:
                index = next(auto_increse)
                kmer_loc[cur_kmer] = index
                cur_loc = index
                # loc_kmer[start]=cur_kmer
                x.append(cur_kmer_feature)
            # kmer muti
            # pdb.set_trace()
            if pre_loc != None:
                u.append(pre_loc)
                v.append(cur_loc)
            pre_kmer = cur_kmer
        if tax_id == 0 or (not tax_id.isdigit()):
            continue
        yield x, (u, v, weight), int(tax_id), {v: k for k, v in kmer_loc.items()}
        u.clear()
        v.clear()
        weight.clear()
        x.clear()
        kmer_loc.clear()


class Y_Handle():
    def __init__(self, file_name):
        self.file_name = file_name

    def muti_label_mode(self):
        self.tree = self._get_ys_tree()
        self.ys_dic, self.level_tax = self._get_one_hot_dic()
        self.tree_prefix_num = len(self.tree.lineage)
        self.max_y = self.tree.get_farthest_leaf()[1]

    def _get_ys_tree(self):
        file_name = self.file_name
        read_generator = ReadGenerator(file_name, "reads")
        res = []
        for read in read_generator.read_Generator():
            tax_id = read.id.split("|kraken:taxid|")[-1].split("_")[0]
            res.append(tax_id)
        if res:
            return ncbi.get_topology(res,  intermediate_nodes=True)
        else:
            raise Exception(f"{file_name} length is zero.")

    def get_genus():
        taxid_levelup_taxid = {}

        def foo(tax_id):
            target = taxid_levelup_taxid.get(tax_id, None)
            if target != None:
                return target
            loc = -1
            temp = ncbi.get_lineage(tax_id)
            target = temp[loc]
            while ncbi.get_rank([target])[target] != "genus":
                loc = loc-1
                target = temp[loc]
            taxid_levelup_taxid[tax_id] = target
            return target
        return foo

    def one_hot_index(self):
        y_dic = {-1: 1}

        def foo(y):
            y_index = y_dic.get(y)
            if y_index == None:
                y_dic[y] = y_dic[-1]
                y_index = y_dic[-1]
                y_dic[-1] = y_dic[-1]+1
            return y_index
        return foo, y_dic

    def _get_one_hot_dic(self):
        _dic = {}
        def dfs(childrens):
            if not childrens:
                return
            one_hot_index_inter, y_dic = self.one_hot_index()
            for tax in childrens:
                one_hot_index_inter(tax.taxid)
                dfs(tax.get_children())
            _dic.update(y_dic)
        childrens = self.tree.get_children()
        dfs(childrens)
        queue = [self.tree]
        res = {"root": self.tree.taxid}
        while queue:
            temp = queue.pop(0)
            level = []
            for tax in temp.get_children():
                level.append(tax.taxid)
                queue.append(tax)
            if level:
                res[temp.taxid]= level
        return _dic,res

    def get_linear_one_hot(self, y):
        _dic = self.ys_dic
        tree_prefix_num = self.tree_prefix_num
        max_y = int(self.max_y)
        temp = ncbi.get_lineage(int(y))[tree_prefix_num:]
        temp = [_dic[taxid] for taxid in temp]
        res = np.pad(temp, (0, int(max_y)-len(temp)),
                     'constant', constant_values=(0))
        return res
