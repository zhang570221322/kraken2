import pickle
import pdb
import matplotlib.pyplot as plt
import arrow
import numpy as np
import sys

plt.rcParams['font.size'] = 20
file_name = sys.argv[1]
lines = open(file_name, "r").readlines()


class STAT:
    def __init__(self) -> None:
        pass


class GENERAL_DATA:
    def __init__(self, number=None, proportion=None) -> None:
        self.number = number
        self.proportion = proportion


class GENUS_LEVEL_DATA:
    def __init__(self, precision, recall, f_measure) -> None:
        self.precision = precision
        self.recall = recall
        self.f_measure = f_measure


def func1(split_str, line):
    temp = line.strip().split(split_str)
    # pdb.set_trace()
    arg1, arg2 = temp[0].strip(), temp[1].strip()[1:-2]
    return float(arg1), float(arg2)


def func2(split_str, line):
    temp = line.split(split_str)[1].strip()[:-2]
    return float(temp)


def func3(line, split_str="processed in"):
    temp = line.split(split_str)[1].split("(")[0].strip()[:-1]
    return float(temp)


def func4(line):
    temp = [temp.split(":") for temp in line[10:].replace(
        "no rank", "no_rank").strip().split(" ")]
    res = {}
    for t in temp:
        if len(t) == 2:
            res[t[0]] = int(t[1])
    return res


def func5(a, b):
    _dic = {}
    for k, _ in a.items():
        _dic[k] = 0
    for k, _ in b.items():
        _dic[k] = 0
    for k, _ in _dic.items():
        if k not in a:
            a[k] = 0
        if k not in b:
            b[k] = 0
    return a, b


def create_stat(lines):
    sequences_classified = GENERAL_DATA(
        *func1("sequences classified", lines[0]))
    missing_Taxonomy = GENERAL_DATA(*func1(
        "sequences taxo is missing in current Taxonomy", lines[1]))
    not_properly_classified = GENERAL_DATA(*func1(
        "sequences is not properly classified", lines[2]))
    classified_genus = GENERAL_DATA(*func1(
        "sequences is properly classified under the all level", lines[3]))
    sequences_unclassified = GENERAL_DATA(*func1(
        "sequences unclassified", lines[4]))
    return sequences_classified, missing_Taxonomy, not_properly_classified, classified_genus, sequences_unclassified


def create_GENUS_LEVEL_DATA(lines, level):
    precision = func2(f"precision at {level} level:", lines[0])
    recall = func2(f"recall at {level} level:", lines[1])
    f_measure = func2(f"f-measure at {level} level:", lines[2])
    return GENUS_LEVEL_DATA(precision, recall, f_measure)


# 总体
kraken2_origin_general = create_stat(lines[7:12])
kraken2_weight_general = create_stat(lines[39:45])
# 冲突次数
conflict_kmers = int(lines[30].split("(kmers:")[1].strip()[:-1])
conflicts_RTL_times = int(lines[31].split("Find conflicts")[
    1].split("times")[:-1][0].strip())
# 时间
kraken2_origin_time = func3(lines[6])
kraken2_weight_conflict_time = func3(lines[31])
kraken2_weight_time = func3(lines[38])

# 属
kraken2_origin_stic_GENUS = create_GENUS_LEVEL_DATA(lines[15:18], "genus")
kraken2_weight_stic_GENUS = create_GENUS_LEVEL_DATA(lines[47:50], "genus")
# 种
kraken2_origin_stic_SPECIES = create_GENUS_LEVEL_DATA(lines[21:24], "species")
kraken2_weight_stic_SPECIES = create_GENUS_LEVEL_DATA(lines[53:56], "species")

# ranks
origin_All_ranks = func4(lines[24])
weight_All_ranks = func4(lines[56])
origin_All_ranks, weight_All_ranks = func5(origin_All_ranks, weight_All_ranks)
data = [kraken2_origin_general, kraken2_weight_general, conflict_kmers,
        conflicts_RTL_times, kraken2_origin_time, kraken2_weight_conflict_time,
        kraken2_weight_time, kraken2_origin_stic_GENUS, kraken2_weight_stic_GENUS,
        kraken2_origin_stic_SPECIES, kraken2_weight_stic_SPECIES, origin_All_ranks,
        weight_All_ranks]
with open(f"{file_name}.data", "wb") as f:
    pickle.dump(data, f)


def bar_number_h(category):
    for rect in category:
        w = rect.get_width()
        plt.text(w, rect.get_y()+rect.get_height() /
                 2, w, ha='left', va='center')


def bar(title, y1, y2, xticks, ylabel, label=['origin', 'weight']):
    # 绘图
    Y_axis = np.arange(len(xticks))

    # plt.barh(Y_axis - 0.2, year2017, 0.4, label = '2017')
    # plt.barh(Y_axis + 0.2, year2018, 0.4, label = '2018')
    if y1:
        b2017 = plt.barh(Y_axis - 0.2, y1, 0.4, label=label[0])
        bar_number_h(b2017)
    if y2:
        b2018 = plt.barh(Y_axis + 0.2, y2, 0.4, label=label[1])
        bar_number_h(b2018)

    plt.yticks(Y_axis, xticks)
    # plt.ylabel("type")
    plt.xlabel(ylabel)
    plt.title(title)
    plt.legend(loc=9)


def my_plot():
    plt.figure(figsize=(57, 30))
    # 总体
    plt.subplot(231)
    y1 = [temp.proportion for temp in kraken2_origin_general]
    y2 = [temp.proportion for temp in kraken2_weight_general]
    ylabel = '(%)'
    xticks = ["sequences_classified", "missing_Taxonomy",
              "not_properly_classified", "classified_allLevel", "sequences_unclassified"]
    bar("kraken2_stic", y1, y2, xticks, ylabel)
    # 冲突次数
    plt.subplot(232)
    y1 = [conflict_kmers, conflicts_RTL_times]
    xticks = ["conflict_kmers", "conflicts_RTL_times"]
    bar("conflict", [], y1, xticks, ylabel, [None, "weight"])

    # 时间
    plt.subplot(233)
    y1 = [kraken2_origin_time, kraken2_weight_conflict_time, kraken2_weight_time]
    y1.append(round(y1[1]+y1[2], 3))
    xticks = ["origin_time", "weight_conflict_time",
              "weight_time", "weight_total"]
    bar("", y1, [], xticks, ylabel, ["time", None])
    # GENUS
    plt.subplot(234)
    y1 = [kraken2_origin_stic_GENUS.precision,
          kraken2_origin_stic_GENUS.recall, kraken2_origin_stic_GENUS.f_measure]
    y2 = [kraken2_weight_stic_GENUS.precision,
          kraken2_weight_stic_GENUS.recall, kraken2_weight_stic_GENUS.f_measure]
    xticks = ["precision", "recall", "f_measure"]
    bar("genus", y1, y2, xticks, ylabel)
    # SPECIES
    plt.subplot(235)
    y1 = [kraken2_origin_stic_SPECIES.precision,
          kraken2_origin_stic_SPECIES.recall, kraken2_origin_stic_SPECIES.f_measure]
    y2 = [kraken2_weight_stic_SPECIES.precision,
          kraken2_weight_stic_SPECIES.recall, kraken2_weight_stic_SPECIES.f_measure]
    xticks = ["precision", "recall", "f-measure"]
    bar("species", y1, y2, xticks, ylabel)
    # all_ranks
    plt.subplot(236)
    xticks = list(weight_All_ranks.keys()) + list((origin_All_ranks.keys()))
    xticks = list(set(xticks))
    y1 = [origin_All_ranks[x] for x in xticks]
    y2 = [weight_All_ranks[x] for x in xticks]
    bar("ranks", y1, y2, xticks, ylabel)
    # time_str = arrow.now().format("YYYY_MM_DD_HH_mm_ss")
    plt.savefig(f"{file_name}.jpg")


my_plot()
