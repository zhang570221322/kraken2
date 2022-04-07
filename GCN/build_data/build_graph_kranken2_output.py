
import pdb
from ete3 import NCBITaxa
ncbi = NCBITaxa()

#%%
def get_total_tqdm(file):
    count=0
    for _ in open(file,"r"):
        count+=1
    return count
def level_traversal(seq_id, tree, out_dic, real,seq_index):
    """层次遍历树
    生成权重x,
    生成邻接矩阵,
    生成y标签,
    """
    temp = [tree]
    weight = []
    taxid_index = {}
    u,v = [],[]
    a = 0
    b = 1
    # 用于Y标签生成
    real_node = None
    try:
        y_list=[]
        while temp:
            node = temp.pop(0)
            if str(node.taxid) == real:
                real_node = node
            weight.append(int(out_dic.get(node.taxid, 0)))
            # weight.append(f"{node.taxid}")
            temp = temp+node.children
            for _ in range(len(node.children)):
                u.append(a)
                v.append(b)
                b = b+1
            y_list.append(node.taxid)
            taxid_index[node.taxid] = a
            a = a+1
        
        Y = taxid_index[real_node.taxid]
        # score = 20
        # Y = [0 for _ in range(len(weight))]
        # score_index = 1
        # if real_node:
        #     temp = real_node
        #     while True:
        #         if score <= 0:
        #             break
        #         index = taxid_index[temp.taxid]
        #         Y[index] = score
        #         score_index += 1
        #         score -= score_index**(1.5)
        #         if temp.is_root():
        #             break
        #         temp = temp.up
        #     _sum = sum(Y)
        #     for i in range(len(Y)):
        #         Y[i] = Y[i]/_sum
        return weight, u,v, Y,y_list,seq_index
    except Exception as e:
        print(f"echo '{str(e)}>>>>>>>>>>{seq_id}' >> error")
    return None


def handle_data(lines_pool, ncbi):
    """
    流式输出: (weight,adjacent_matrix,Y)
    """
    for seq_index,line in enumerate(lines_pool):
        line = line.strip()
        line_split = line.split('\t')
        seq_id = line_split[1]
        real = seq_id.split("|")[-1].split("_")[0]
        # predi = line_split[2]
        out_dic = {}
        for i in line_split[4].strip().split(' '):
            _i = i.split(":")
            value = int(_i[1])
            # 如果为非数字
            if not _i[0].isdigit():
                continue
            key = int(ncbi._translate_merged([str(_i[0])])[0].pop())
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
                               real,seq_index)
        if not temp:
            continue
        yield temp
def get_feature(file):
    f = open(file,"r")
    for i in handle_data(f,ncbi):
        yield i
    f.close()

#%%
if __name__=="__main__":
    lines =open("/home/yegpu/zwl/kraken2/GCN/raw/0.fastq.origin.output","r").readlines()
    for i in handle_data(lines,ncbi):
        print(i)
        break
# %%
