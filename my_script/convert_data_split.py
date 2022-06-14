import numpy as np
from ete3 import NCBITaxa
import random
import multiprocessing
import pdb
print("nohup python3 convert_data.py > log 2>&1 &")
# 公用临时变量


def choice(probability):
    """
    P(范围值=True)=probability
    """
    if random.random() <= probability:
        return True


# def get_merged_taxid(taxid, merged_dic, ncbi):
#     """
#     返回更改后的taxid
#     taxid:int
#     """
#     if taxid in merged_dic:
#         return merged_dic[taxid]
#     merged_dic.update(ncbi._translate_merged([str(taxid)])[1])
#     return merged_dic[taxid]


def level_traversal(seq_id, tree, out_dic, real, out_order_list, merged_dic, ncbi):
    """层次遍历树
    生成权重x,
    生成邻接矩阵,
    生成y标签,
    """
    default_N = 220
    temp = [tree]
    weight = []
    taxid_index = {}
    # 用于创建邻接矩阵
    # adjacent_matrix = [[0 for _ in range(default_N)] for _ in range(default_N)]
    adjacent_matrix = []
    a = 0
    b = 1
    # 用于Y标签生成
    real_node = None
    try:
        while temp:
            node = temp.pop(0)
            if str(node.taxid) == real:
                real_node = node
            weight.append(int(out_dic[node.taxid]))
            # weight.append(f"{node.taxid}")
            temp = temp+node.children
            for _ in range(len(node.children)):
                adjacent_matrix.append((a, b))
                b = b+1
            taxid_index[node.taxid] = a
            a = a+1
        # 补全weight
        # for _ in range(default_N-len(weight)):
        #     weight.append(0)
        score = 20
        Y = [0 for _ in range(len(weight))]
        score_index = 1
        if real_node:
            temp = real_node
            while True:
                if score <= 0:
                    break
                index = taxid_index[temp.taxid]
                Y[index] = score
                score_index += 1
                score -= score_index**(1.5)
                if temp.is_root():
                    break
                temp = temp.up
            _sum = sum(Y)
            for i in range(len(Y)):
                Y[i] = Y[i]/_sum
            for i in range(len(out_order_list)):
                taxid = out_order_list[i]
                if taxid == 0 or taxid == -2:
                    continue
                order_id = taxid_index[taxid] + 1
                out_order_list[i] = order_id
            return weight, out_order_list, adjacent_matrix, Y
    except Exception as e:
        print(str(e))
        print(seq_id)
        print(">>>>>>>>>")
        raise e
    return None


def handle_data(lines_pool, merged_dic, ncbi):
    """
    流式输出: (weight,adjacent_matrix,Y)
    """
    for line in lines_pool:
        line = line.strip()
        line_split = line.split('\t')
        line_type = int(line_split[0])
        #  排除dberror和error以及未分类的
        if line_type in [3, 5, 6]:
            continue
        #  只选择一部分经典标签
        # if line_type in [0, 1]:
        #     if choice(0.8):
        #         continue
        seq_id = line_split[2]
        real = seq_id.split("|")[-1].split("_")[0]
        # predi = line_split[3]
        out_dic = {}
        out_order_list = []
        for i in line_split[5].strip().split(' '):
            _i = i.split(":")
            value = int(_i[1])
            # 如果为非数字
            if not _i[0].isdigit():
                [out_order_list.append(-2) for _ in range(value)]
                continue
            key = int(ncbi._translate_merged([str(_i[0])])[0].pop())
            [out_order_list.append(key) for _ in range(value)]
            if key == 0:
                continue
            value_count = out_dic.get(key, None)
            if value_count != None:
                out_dic[key] = out_dic[key]+value
            else:
                out_dic[key] = value
        out_list = list(out_dic.keys())
        # 排除过少的命中组合
        if len(out_list) <= 2:
            continue
        # 排除一个也没命中的
        if int(real) not in out_list:
            continue
        tree = ncbi.get_topology(out_list, intermediate_nodes=False)
        temp = level_traversal(seq_id, tree, out_dic,
                               real, out_order_list, merged_dic, ncbi)
        print(temp)
        if not temp:
            continue
        yield temp


# 公用输出变量

def data_iter(batch, data):
    """
    # 均分data
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]


def func(index, lines_pool):
    x1 = []
    x2 = []
    adjacent_matrixs = []
    Ys = []
    merged_dic = {}
    ncbi = NCBITaxa()
    for weight, out_order_list, adjacent_matrix, Y in handle_data(lines_pool, merged_dic, ncbi):
        x1.append(weight)
        x2.append(out_order_list)
        adjacent_matrixs.append(adjacent_matrix)
        Ys.append(Y)
    # np.save(f"x1_{index}", x1)
    # np.save(f"x2_{index}", x2)
    # np.save(f"adjacent_matrixs_{index}", adjacent_matrixs)
    # np.save(f"Y_{index}", Ys)


result = "/data/home/wlzhang/classfication/kraken2_db/2021_11_01/library/stimulation/test/1.output.typed"
lines = open(
    result, "r").readlines()
processing_nums = 1

for index, lines_pool in enumerate(data_iter(processing_nums, lines)):
    print(f"exec process {index} len(lines_pool)={len(lines_pool)}")
    func(index, lines_pool)


def get_merged_taxid(taxid, merged_dic, ncbi):
    """
    返回更改后的taxid
    taxid:int
    """
    if taxid in merged_dic:
        return merged_dic[taxid]
    merged_dic.update(ncbi._translate_merged([str(taxid)])[1])
    return merged_dic[taxid]


def get_taxid_index(taxid_index, taxid, merged_dic, ncbi):
    """
    {}.get(a,func) , 会在get前执行func
    """
    temp = taxid_index.get(taxid, None)
    if temp != None:
        return temp
    else:
        return get_merged_taxid(taxid, merged_dic, ncbi)
