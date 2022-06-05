
import pdb
from common import log_time, AUTO_INCREMENT
from build_data.graph_model import   ReadGenerator,  feature_str,get_taxid_kraken2
import copy
from ete3 import NCBITaxa
import numpy as np
from build_data.load_ncbi_taxinfo import get_id_path,taxonomy
 

ncbi = NCBITaxa()

 

def data_iter(batch, data):
    """
    # 均分data
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]



def get_feature(file_name):
    read_generator = ReadGenerator(file_name, "reads")
    _kmer_loc={}
    x,u, v, weight = [],[], [], []
 
    for read in read_generator.read_Generator():
        tax_id = 0

        tax_id = get_taxid_kraken2(read.id)
        for kmer, kmer_read_index, _ in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            temp = _kmer_loc.get(kmer,None)
            if temp:
                temp.append(start)
            else:
                _kmer_loc[kmer]=[start]
        for kmer, kmer_read_index, cur_kmer_feature in read_generator.Kmer_index_Generator(read):
            x.append(cur_kmer_feature)
            u_index, _ = kmer_read_index
            v_list = _kmer_loc[kmer]
            for v_index in v_list:
                u.append(u_index)
                v.append(v_index)
                weight.append(abs(u_index-v_index))
        if tax_id == 0 or (not tax_id.isdigit()):
            continue
        yield x, (u, v, weight), int(tax_id) 
        u.clear()
        v.clear()
        weight.clear()
        _kmer_loc.clear()
        x.clear()
def get_feature2(file_name):
    read_generator = ReadGenerator(file_name, "reads")
    x, u, v, weight = [], [], [], []
    kmer_loc = {}
    # loc_kmer={}
    for read in read_generator.read_Generator():
        tax_id = 0
        pre_kmer = ""
        auto_increse = AUTO_INCREMENT()
        tax_id = get_taxid_kraken2(read.id)
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
        
def get_feature3(file_name):
    read_generator = ReadGenerator(file_name, "fasta2fastq")
    x, u, v, weight = [], [], [], []
    kmer_loc = {}
    # loc_kmer={}
    for read in read_generator.read_Generator():
        tax_id = 0
        pre_kmer = ""
        auto_increse = AUTO_INCREMENT()
        tax_id = get_taxid_kraken2(read.id)
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
    def __init__(self, file_name,type="reads"):
        self.file_name = file_name
        self.count=0
        self.type=type

    def muti_label_mode(self,mode=1):
        self.tree = self._get_ys_tree()
        if mode == 1:
            self.ys_dic, self.level_tax,self.level_dim = self._get_one_hot_dic()
        elif mode == 2:
            self.ys_dic, self.level_tax,self.level_dim = self._get_one_hot_dic2()
        self.max_y = 7

    def _get_ys_tree(self):
        file_name = self.file_name
        read_generator = ReadGenerator(file_name, self.type)
        res = set()
        all_id_path={}
        for read in read_generator.read_Generator():
            tax_id = get_taxid_kraken2(read.id)
            id_paths =get_id_path(tax_id,*taxonomy)
            for id_path in id_paths:
                res.add(id_path)
            if int(tax_id) not in all_id_path:
                all_id_path[int(tax_id)]=(id_paths)
        if res:
            self.all_id_path = all_id_path
            return ncbi.get_topology(res,  intermediate_nodes=False)
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
    def _get_one_hot_dic2(self):
        res = {0:set(),1:set(),2:set(),3:set(),4:set(),5:set(),6:set()}

        for id_path in self.all_id_path.values():
            for i in range(1,len(id_path)):
                level = i-1
                res[level].add(id_path[i])
        _dim={}
        _dic={}
        _level_tax={}
        for level,v in res.items():
            _dim[level]= len(v)+1
            one_hot_index_inter, y_dic = self.one_hot_index()
            for tax in v:
                index = one_hot_index_inter(tax)
                _level_tax[f"{level}_{index}"]= tax
            _level_tax[f"{level}_0"] = 0
            _dic.update(y_dic)
        
        return _dic,_level_tax,_dim

    def _get_one_hot_dic(self):
        _dic = {}
        _dim = {}
        def dfs(childrens,level):
            if not childrens:
                return
            one_hot_index_inter, y_dic = self.one_hot_index()
            for tax in childrens:
                one_hot_index_inter(tax.taxid)
                dfs(tax.get_children(),level+1)
            _dic.update(y_dic)
            if level in _dim:
                _dim[level] = max(_dim[level],y_dic[-1])
            else:
                _dim[level] = y_dic[-1]
        childrens = self.tree.get_children()
        dfs(childrens,0)
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
        return _dic,res,_dim

    def get_linear_one_hot(self, y):
        _dic = self.ys_dic
        max_y=self.max_y
        id_path = self.all_id_path[int(y)][1:]
        # try:
        
        # self.count+=1
        # print(self.count)
        # if self.count==180:
        #     pdb.set_trace()            
        temp = []
        for taxid in id_path:
            # if taxid ==0:
            #     temp.append(0)
            # else:
            temp.append(_dic[taxid])
        # except KeyError as err:
        #     pdb.set_trace()
        #     print(id_path)
        #     print(y)
        #     raise err

        res = np.pad(temp, (0, int(max_y)-len(temp)),
                     'constant', constant_values=(0))
        return res,id_path
