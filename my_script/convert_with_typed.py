
"""统计kraken2的结果
"""
# %%

from enum import Enum
from ete3 import NCBITaxa
from collections import ChainMap
import os
import argparse

# %%
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


class Type(Enum):
    total_assegned_s = 0
    total_assegned_g = 1
    other = 2
    total_unclassified = 3
    need_improve = 4
    total_dberror = 5
    total_error = 6


class stats:
    total_sequences = 0
    total_classified = 0
    total_assegned_g = 0  # TP
    total_assegned_s = 0
    total_unclassified = 0  # FN
    total_dberror = 0
    total_error = 0
    need_improve = 0
    other = 0

    def __str__(self) -> str:
        return f"total_sequences: {self.total_sequences}\n\
            total_classified: {self.total_classified}\n\
            total_assegned_g:{self.total_assegned_g}\n\
            total_assegned_s:{self.total_assegned_s} \n\
            total_unclassified:{self.total_unclassified}\n\
            total_dberror:{self.total_dberror}\n\
            total_error:{self.total_error }\n\
            need_improve:{self.need_improve}\n\
            other:{self.other}\n"


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
    try:
        t1 = temp[0]
        t2 = temp[1]
    except IndexError:
        return -1
    x1 = x2 = 0
    i2 = 0
    if len(t1) > len(t2):
        x1, x2 = len(t1), len(t2)
        t1, t2 = temp[0], temp[1]
    else:
        x1, x2 = len(t2), len(t1)
        t1, t2 = temp[1], temp[0]
    for i1 in range(x1):
        if i2 >= x2:
            return t2[-1]
        if t1[i1] != t2[i2]:
            return t1[i2]
        i2 += 1
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


# %%
#  原始数据
right = combined["mapping_file"]
result = combined["output_file"]
typed = open(f"{result}.typed", "w")
typed_result = open(f"{result}.typed.result", "w")


# 正确的序列id和对应的真是taxid.
data_result = {}
for line in open(right, 'r'):
    title = line.strip()
    split_title = title.split("\t")
    ID = split_title[0]
    taxid = split_title[1]
    data_result[ID] = taxid
print("加载参考数据完成")


def compute1(result, stats, data_result):
    with open(result, "r") as f:
        line = f.readline().strip()
        while line:
            try:
                line_split = line.split('\t')
                is_classification = line_split[0] == 'C'
                seq_id = line_split[1]
                predi = line_split[2]
                out_list = set([i.split(':')[0]
                                for i in line_split[4].split(' ')])
                stats.total_sequences += 1
                if not is_classification:  # 如果一个也没命中
                    stats.total_unclassified += 1
                    typed.write(f"{Type.total_unclassified.value}\t{line}\n")
                    line = f.readline().strip()
                    continue
                if int(predi) == int(data_result[seq_id]):  # 如果预测命中
                    rank = ncbi.get_rank([predi])[int(predi)]
                    if rank in species_types:
                        stats.total_assegned_g += 1
                        stats.total_assegned_s += 1
                        typed.write(f"{Type.total_assegned_s.value}\t{line}\n")
                    elif rank in genus_types:
                        stats.total_assegned_g += 1
                        typed.write(f"{Type.total_assegned_g.value}\t{line}\n")
                    elif real in out_list:
                        stats.need_improve += 1
                        typed.write(f"{Type.total_assegned_g.value}\t{line}\n")
                    else:
                        stats.other += 1
                        typed.write(f"{Type.other.value}\t{line}\n")
                else:  # 如果不对, 看是否是上升到LCA了.
                    real = data_result[seq_id].strip()
                    common_node = ncbi.get_common_ancestor(
                        real, predi)
                    if common_node:
                        if str(common_node) == predi:
                            rank = ncbi.get_rank([predi])[int(predi)]
                            if rank in species_types:
                                stats.total_assegned_g += 1
                                stats.total_assegned_s += 1
                                typed.write(
                                    f"{Type.total_assegned_s.value}\t{line}\n")
                            elif rank in genus_types:
                                stats.total_assegned_g += 1
                                typed.write(
                                    f"{Type.total_assegned_g.value}\t{line}\n")
                            elif real in out_list:
                                stats.need_improve += 1
                                typed.write(
                                    f"{Type.total_assegned_g.value}\t{line}\n")
                            else:
                                stats.other += 1
                                typed.write(f"{Type.other.value}\t{line}\n")
                        # 如果真实taxid存在kmer的命中列表,但是没有被选中.
                        elif real in out_list:
                            stats.need_improve += 1
                            typed.write(
                                f"{Type.total_assegned_g.value}\t{line}\n")
                        else:
                            stats.total_error += 1
                            typed.write(
                                f"{Type.total_error.value}\t{line}\n")
                    else:
                        stats.total_error += 1
                        typed.write(
                            f"{Type.total_error.value}\t{line}\n")
                line = f.readline().strip()
            except KeyError:
                stats.total_dberror += 1
                typed.write(
                    f"{Type.total_dberror.value}\t{line}\n")
                line = f.readline().strip()


print("start")
compute1(result, stats, data_result)
typed.close()
stats.total_classified = stats.total_sequences - stats.total_unclassified
precision_g = 1.0 * stats.total_assegned_g / stats.total_classified
precision_s = 1.0 * stats.total_assegned_s / stats.total_classified
recall_g = 1.0 * stats.total_assegned_g / stats.total_sequences
recall_s = 1.0 * stats.total_assegned_s / stats.total_sequences
f_mesure_g = (2.0 * precision_g * recall_g) / (precision_g +
                                               recall_g) if (recall_g != 0 and precision_g != 0) else -1.0
f_mesure_s = (2.0 * precision_s * recall_s) / (precision_s +
                                               recall_s)if (recall_s != 0 and precision_s != 0)else -1.0
typed_result.write(stats.__str__(stats))
typed_result.write("GENERAL DATA\n")
typed_result.write("  %d sequences classified (%.4f%%)\n" % (
    stats.total_classified,
    stats.total_classified * 100.0 / stats.total_sequences))
typed_result.write("  %d sequences unclassified (%.4f%%)\n\n" % (
    stats.total_unclassified,
    stats.total_unclassified * 100.0 / stats.total_sequences))
typed_result.write("GENUS LEVEL DATA\n")
typed_result.write("  %d sequences classified at genus level.\n" % (
    stats.total_assegned_g))
typed_result.write("  precision at genus level: %.4f%%.\n" %
                   (100.0 * precision_g))
typed_result.write("  recall at genus level: %.4f%%.\n" % (100.0 * recall_g))
typed_result.write("  f-measure at genus level: %.4f%%.\n\n" %
                   (100.0 * f_mesure_g))
typed_result.write("SPECIES LEVEL DATA\n")
typed_result.write("  %d sequences classified at species level.\n" %
                   (stats.total_assegned_s))
typed_result.write("  precision at species level: %.4f%%.\n" %
                   (100.0 * precision_s))
typed_result.write("  recall at species level: %.4f%%.\n" % (100.0 * recall_s))
typed_result.write("  f-measure at species level: %.4f%%.\n" %
                   (100.0 * f_mesure_s))
typed_result.close()
# %%
