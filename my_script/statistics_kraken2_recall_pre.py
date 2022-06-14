
"""统计kraken2的结果
"""
# %%
import pdb
from ete3 import NCBITaxa
from collections import ChainMap
import os
import argparse

# %%
print("Test")
defaults = {
    "mapping_file": "Empty",
    "output_file": "Empty"
}
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mapping_file')
parser.add_argument('-o', '--output_file')
namespace = parser.parse_args()
command_line_args = {k: v for k, v in vars(namespace).items() if v}
combined = ChainMap(command_line_args, os.environ, defaults)

ncbi = NCBITaxa()

species_types = [
    "species",
    "subspecies",
    "variety",
    "subvariety",
    "form",
    "subform",
    "strain",
    "no rank",
    "species subgroup",
    "species group",
    "forma specialis"]
genus_types = [
    "genus",
    "subgenus",
    "section",
    "subsection",
    "series",
    "subseries"]
others_types = [
    "root",
    "superkingdom",
    "kingdom",
    "subkingdom ",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "superorder",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "biotype"]

print("加载完成数据库")


def a(tax1: str, tax2: str):
    """返回tax1和tax2的最低共同祖先
    """
    if tax1 == tax2:
        return int(tax1)
    temp = list(ncbi.get_lineage_translator([tax1, tax2]).values())
    common = set(temp[0]) & set(temp[1])
    if common:
        return list(common)[-1]
    else:
        return None


ncbi.get_common_ancestor = a


def IsSpecies(rank: str) -> bool:
    is_species = False
    if rank in species_types:
        is_species = True
    return is_species


def IsGenus(rank):
    is_genus = False
    if rank in genus_types:
        is_genus = True
    return is_genus


def IsOther(rank):
    is_other = False
    if rank in others_types:
        is_other = True
    return is_other


def compute1(result, stats, data_result):
    with open(result, "r") as f:
        line = f.readline()
        while line:
            try:
                line_split = line[:-1].split('\t')
                # pdb.set_trace()
                is_classification = line_split[0] == 'C'
                seq_id = line_split[1]
                predi = line_split[2]
                out_list = set([i.split(':')[0]
                                for i in line_split[4].split(' ')])
                line = f.readline()
                stats.total_sequences += 1
                # 如果一个也没命中
                if not is_classification:
                    stats.total_unclassified += 1
                    continue
                # 如果预测命中
                if int(predi) == int(data_result[seq_id]):
                    print(seq_id)
                    rank = ncbi.get_rank([predi])[int(predi)]
                    if rank in species_types:
                        stats.total_assegned_g += 1
                        stats.total_assegned_s += 1
                    elif rank in genus_types:
                        stats.total_assegned_g += 1
                # 如果不对, 看是否是上升到LCA了.
                else:
                    real = data_result[seq_id]
                    common_node = ncbi.get_common_ancestor(
                        real, predi)
                    if common_node:
                        if str(common_node) == predi:
                            print(seq_id)
                            rank = ncbi.get_rank([predi])[int(predi)]
                            if rank in species_types:
                                stats.total_assegned_g += 1
                                stats.total_assegned_s += 1
                            elif rank in genus_types:
                                stats.total_assegned_g += 1
                        elif common_node == -1:
                            stats.total_dberror += 1
                    else:
                        stats.total_error += 1
                    # 如果真实taxid存在kmer的命中列表,但是没有被选中.
                    if data_result[seq_id] in out_list:
                        stats.need_improve += 1
            except KeyError:
                stats.total_dberror += 1


def compute2(result, stats, data_result):
    with open(result, "r") as f:
        line = f.readline()
        while line:
            try:
                line_split = line[:-1].split('\t')
                is_classification = line_split[0] == 'C'
                seq_id = line_split[1]
                predi = line_split[2]
                out_list = set([i.split(':')[0]
                                for i in line_split[4].split(' ')])
                real = data_result[seq_id]
                line = f.readline()
                stats.total_sequences += 1

                # 如果一个也没命中
                if not is_classification:
                    stats.total_unclassified += 1
                    continue
                rank = ncbi.get_rank([predi])[int(predi)]
                # 预测的节点只有是真实节点的父亲或等于真实节点. 才算预成功.
                if ncbi.get_topology([int(predi), int(real)]).name == str(predi):
                    # 检查是否精确到种
                    if IsSpecies(rank):
                        stats.total_assegned_g += 1
                        stats.total_assegned_s += 1
                        continue
                    # 如果精确到属
                    if IsGenus(rank):
                        stats.total_assegned_g += 1
                        continue
                else:
                    stats.total_error += 1
                # 如果真实taxid存在kmer的命中列表,但是没有被选中.
                # 精确选中
                if data_result[seq_id] in out_list:
                    stats.need_improve += 1
            except KeyError:
                stats.total_dberror += 1


# %%
#  原始数据
right = combined["mapping_file"]
result = combined["output_file"]
# 正确的序列id和对应的真是taxid.
data_result = {}
for line in open(right, 'r'):
    title = line.strip()
    split_title = title.split("\t")
    ID = split_title[0]
    taxid = split_title[1]
    data_result[ID] = taxid
print("加载参考数据完成")


class stats:
    total_sequences = 0
    total_classified = 0
    total_assegned_g = 0  # TP
    total_assegned_s = 0
    total_unclassified = 0  # FN
    total_dberror = 0
    total_error = 0
    need_improve = 0


print("start")
compute1(result, stats, data_result)
stats.total_classified = stats.total_sequences - stats.total_unclassified
precision_g = 1.0 * stats.total_assegned_g / stats.total_classified
precision_s = 1.0 * stats.total_assegned_s / stats.total_classified
recall_g = 1.0 * stats.total_assegned_g / stats.total_sequences
recall_s = 1.0 * stats.total_assegned_s / stats.total_sequences
f_mesure_g = (2.0 * precision_g * recall_g) / (precision_g +
                                               recall_g) if (recall_g != 0 and precision_g != 0) else -1.0
f_mesure_s = (2.0 * precision_s * recall_s) / (precision_s +
                                               recall_s)if (recall_s != 0 and precision_s != 0)else -1.0

print("GENERAL DATA")
print("  %d sequences classified (%.4f%%)" % (
      stats.total_classified,
      stats.total_classified * 100.0 / stats.total_sequences))
print("  %d sequences unclassified (%.4f%%)" % (
      stats.total_unclassified,
      stats.total_unclassified * 100.0 / stats.total_sequences))
print("GENUS LEVEL DATA")
print("  %d sequences classified at genus level." % (
      stats.total_assegned_g))
print("  precision at genus level: %.4f%%." % (100.0 * precision_g))
print("  recall at genus level: %.4f%%." % (100.0 * recall_g))
print("  f-measure at genus level: %.4f%%." % (100.0 * f_mesure_g))
print("SPECIES LEVEL DATA")
print("  %d sequences classified at species level." %
      (stats.total_assegned_s))
print("  precision at species level: %.4f%%." % (100.0 * precision_s))
print("  recall at species level: %.4f%%." % (100.0 * recall_s))
print("  f-measure at species level: %.4f%%." % (100.0 * f_mesure_s))

# %%
