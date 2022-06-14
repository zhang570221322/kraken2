import pdb
import pickle
import matplotlib.pyplot as plt
import numpy as np
import sys
import os


plt.rcParams['font.size'] = 20
np_mean = np.mean


def f(data):
    return np.round(np_mean(data, axis=0), 4)


np.mean = f


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


files = os.listdir(sys.argv[1])
my_files = [sys.argv[1]+"/"+file for file in files if file.endswith(".data")]
data = []

for file in my_files:
    data.append(pickle.load(open(file, 'rb')))


def func1(number):
    r = range(len(my_files))
    return [[i.proportion for i in temp] for temp in [data[i][number] for i in r]]


def func2(number):
    r = range(len(my_files))
    return [data[i][number] for i in r]


def func3(number):
    r = range(len(my_files))
    return [[temp.precision, temp.recall, temp.f_measure] for temp in [data[i][number] for i in r]]


def func4(number):
    r = range(len(my_files))
    _dic = {}
    for temp in [data[i][number] for i in r]:
        for k, v in temp.items():
            _dic[k] = 0
    for temp in [data[i][number] for i in r]:
        for k, v in temp.items():
            _dic[k] += v
    for k, v in data[0][number].items():
        _dic[k] = v/len(my_files)
    return _dic


kraken2_origin_general = func1(0)
kraken2_weight_general = func1(1)
# 冲突次数
conflict_kmers = np.array(func2(2), dtype=np.float).mean()
conflicts_RTL_times = np.array(func2(3), dtype=np.float).mean()
# 时间
kraken2_origin_time = np.array(func2(4), dtype=np.float).mean()
kraken2_weight_conflict_time = np.array(func2(5), dtype=np.float).mean()
kraken2_weight_time = np.array(func2(6), dtype=np.float).mean()
# 属
kraken2_origin_stic_GENUS = func3(7)
kraken2_weight_stic_GENUS = func3(8)
# 种
kraken2_origin_stic_SPECIES = func3(9)
kraken2_weight_stic_SPECIES = func3(10)
# ranks
origin_All_ranks = func4(11)
weight_All_ranks = func4(12)


def bar_number_h(category):
    for rect in category:
        w = rect.get_width()
        plt.text(w, rect.get_y()+rect.get_height() /
                 2, w, ha='left', va='center')


def bar(title, y1, y2, xticks, ylabel, label=['origin', 'weight']):
    # 绘图
    Y_axis = np.arange(len(xticks))

    if y1 if isinstance(y1, list) else y1.all():
        b2017 = plt.barh(Y_axis - 0.2, y1, 0.4, label=label[0])
        bar_number_h(b2017)
    if y2 if isinstance(y2, list) else y2.all():
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
    y1 = list(np.mean(kraken2_origin_general))
    y2 = list(np.mean(kraken2_weight_general))

    ylabel = '(%)'
    xticks = ["sequences_classified", "missing_Taxonomy",
              "not_properly_classified", "classified_allLevel", "sequences_unclassified"]
    bar("kraken2_stic", y1, y2, xticks, ylabel)
    # 冲突次数
    plt.subplot(232)
    y1 = np.array([conflict_kmers, conflicts_RTL_times])
    xticks = ["conflict_kmers", "conflicts_RTL_times"]
    bar("conflict", [], y1, xticks, ylabel, [None, "weight"])

    # 时间
    plt.subplot(233)
    y1 = [kraken2_origin_time, kraken2_weight_conflict_time, kraken2_weight_time]
    y1.append(round(y1[1]+y1[2], 3))
    y1 = np.array(y1)
    xticks = ["origin_time", "weight_conflict_time",
              "weight_time", "weight_total"]
    bar("", y1, [], xticks, ylabel, ["time", None])
    # GENUS
    plt.subplot(234)
    y1 = np.mean(kraken2_origin_stic_GENUS)
    y2 = np.mean(kraken2_weight_stic_GENUS)
    xticks = ["precision", "recall", "f_measure"]
    bar("genus", y1, y2, xticks, ylabel)
    # SPECIES
    plt.subplot(235)
    y1 = np.mean(kraken2_origin_stic_SPECIES)
    y2 = np.mean(kraken2_weight_stic_SPECIES)
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
    plt.savefig(f"{sys.argv[1]}/total.jpg")


my_plot()
