import pickle
from ete3 import NCBITaxa
import pdb
import sys

from ete3.coretype.tree import Tree
ncbi = NCBITaxa()
file_name = sys.argv[1]
# filter1 = [line.strip() for line in open(sys.argv[2], "r")]
filter1 = []
tree_print = True
meta_print = True
temp_print = print
a = False

_sum_array = []
_predi_sum_array = []
_real_sum_array = []
_predi_parent_array = []
_real_parent_array = []


def my_get(tree, down_id):
    if down_id == int(tree.name):
        return tree
    for temp in tree.get_descendants():
        if int(temp.name) == down_id:
            return temp


def up(tree, out_list):
    if not tree.up and int(tree.name) in out_list:
        return int(tree.name)
    while tree and tree.up and int(tree.up.name) not in out_list:
        tree = tree.up
    if tree and tree.up:
        return int(tree.up.name)
    else:
        return None


def print(t, p):
    if(t):
        temp_print(p)


for line in open(file_name, "r"):
    line_split = line.strip().split('\t')
    if line_split[1] in filter1:
        continue
    is_classification = line_split[0] == 'C'
    if not is_classification:
        continue

    seq_id = line_split[1]
    predi = int(line_split[2])
    right_tax = int(seq_id.split('|')[-1].split('_')[0])
    out_dic = {}
    for i in line_split[4].strip().split(' '):
        _i = i.split(":")
        if not _i[0].isdigit():
            continue
        key = int(_i[0])
        value = int(_i[1])
        if key == 0:
            continue
        temp = out_dic.get(key)
        if temp:
            out_dic[key] = out_dic[key]+value
        else:
            out_dic[key] = value
    # 添加真实的标签
    _total = sum(out_dic.values())
    if right_tax not in out_dic:
        out_dic[right_tax] = -1
    out_list = list(out_dic.keys())

    real_sum_pro = round(out_dic[right_tax]/_total, 3)
    predi_sum_pro = round(out_dic[predi]/_total, 3)
    if len(out_list) == 1:
        real_parent_pro = 0
        predi_parent_pro = 1.0
        # print(meta_print,
        #       f"seq_id: {seq_id}\nsum: {_total},\n真实{right_tax}占总共命中:{real_sum_pro}, 父亲:{real_parent_pro},\n预测{predi}占总共命中:{predi_sum_pro}, 父亲:{predi_parent_pro}")
        # print(meta_print, f"-{out_dic[out_list[0]]}, {out_list[0]}")
        _sum_array.append(_total)
        _predi_sum_array.append(predi_sum_pro)
        _real_sum_array.append(real_sum_pro)
        _predi_parent_array.append(predi_parent_pro)
        _real_parent_array.append(real_parent_pro)
        continue
    tree = ncbi.get_topology(out_list,  intermediate_nodes=False)
    real_parent = up(my_get(tree, right_tax), out_list)
    if not real_parent:
        continue
    real_parent_pro = out_dic[right_tax] / out_dic[real_parent]
    predi_parent = up(my_get(tree, predi), out_list)
    if not predi_parent:
        continue
    predi_parent_pro = out_dic[predi] / out_dic[predi_parent]

    # print(meta_print,
    #   f"seq_id: {seq_id}\nsum: {_total},\n真实{right_tax}占总共命中:{real_sum_pro}, 父亲{real_parent}:{real_parent_pro},\n预测{predi}占总共命中:{predi_sum_pro}, 父亲{predi_parent}:{predi_parent_pro}")
    _sum_array.append(_total)
    _predi_sum_array.append(predi_sum_pro)
    _real_sum_array.append(real_sum_pro)
    _predi_parent_array.append(predi_parent_pro)
    _real_parent_array.append(real_parent_pro)
    # tax2names = ncbi.annotate_tree(
    #     tree, taxid_attr="name", tax2name=out_dic)[0]
    # for k, v in out_dic.items():
    #     tax2names[k] = ","+str(v)
    # nothing = ncbi.annotate_tree(
    #     tree, taxid_attr="name", tax2name=tax2names)
    # print(tree_print,  tree.get_ascii(
    #     attributes=["sci_name", "rank", "taxid"]))
t = [_sum_array, _predi_sum_array, _real_sum_array,
     _predi_parent_array, _real_parent_array]
with open(f"{file_name}.data", "wb") as f:
    pickle.dump(t, f)
#                                                                                                           / 496866: 1, 340099: 1
#                                                                      /68295:1,186814:1,1754:1
#                                                                     |                                     \2636821:2 ,399726:1
# -131567:1, 2:1, 1783272:1, 1239:1, 186801:1
#                                                                     |                                                  / 36827:1, 591968:1, 498213:1
#                                                                     |             / 31979:1,1485:1,1491:1
#                                                                      \ 186802:2                                   \36826:2, 498214:1
#                                                                                  |
#                                                                                   \ 31984:2 , 2831443:1, 35701:1,498761:1
taxs = {8: [['340099'], ['399726'], ['36827', '36826'], ['498761']],
        7: [['496866', '2636821'], ['1491'], ['35701']],
        6: [['1754'], ['1485'], ['2831443']],
        5: [['186814'], ['31979', '31984']],
        10: [['498213']],
        9: [['591968'], ['498214']],
        4: [['68295', '186802']],
        3: [['186801']],
        2: [['1239']],
        1: [['1783272']],
        0: [['2']]}
