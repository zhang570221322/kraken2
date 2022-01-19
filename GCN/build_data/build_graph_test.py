import multiprocessing
import os
from scipy import sparse
import numpy as np
from graph_model import feature_space_init, ReadGenerator, DEFAULT_K
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


def set_value(adj_matrix, feature_index, read_kmers, end, tar_value, feature_space):
    fea_target_feature = feature_space[read_kmers[end]][0]
    value = abs(adj_matrix[feature_index][fea_target_feature])
    if value == 0:
        adj_matrix[feature_index][fea_target_feature] = tar_value
    elif abs(value) > abs(tar_value):
        adj_matrix[feature_index][fea_target_feature] = tar_value
    elif abs(value) == abs(tar_value) and value < 0:
        adj_matrix[feature_index][fea_target_feature] = tar_value


def get_feature(file_name):
    feature_space = copy.deepcopy(feature_space_init)
    adj_matrix = np.zeros((4**DEFAULT_K, 4**DEFAULT_K), dtype=np.int)
    read_generator = ReadGenerator(file_name, "fasta")
    tax_id = 0
    for read in read_generator.read_Generator():
        read_kmers = []
        tax_id = read.id.split("|kraken:taxid|")[-1].split(" ")[0]
        for kmer, _ in read_generator.Kmer_index_Generator(read):
            feature_space[kmer][1] += 1
            read_kmers.append(kmer)
        for kmer, kmer_read_index in read_generator.Kmer_index_Generator(read):
            start, end = kmer_read_index
            feature_index = feature_space[kmer][0]
            for i in range(1, 30):
                if start-i > 0:
                    set_value(adj_matrix, feature_index,
                              read_kmers, start-i, -i, feature_space)
                if end + i <= len(read.seq):
                    set_value(adj_matrix, feature_index,
                              read_kmers, start+i, i, feature_space)
        read_kmers.clear()
    return list(feature_space.values()), sparse.csc_matrix(adj_matrix), tax_id


def func(index, files):
    feature = []
    adjacent_matrixs = []
    Ys = []
    for file_name in files:
        file_dir = f"{prefix}/{file_name}"
        x, adj, y = get_feature(file_dir)
        feature.append(x)
        adjacent_matrixs.append(adj)
        Ys.append(y)
    np.save(f"{our_dir}/x_{index}", feature)
    np.save(f"{our_dir}/adjacent_matrixs_{index}", adjacent_matrixs)
    np.save(f"{our_dir}/Y_{index}", Ys)


our_dir = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/dp2/data/out_data2"
prefix = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/Bacteria_artificial/reference_delete_n"
all_files = os.listdir(prefix)
all_files = sorted(all_files)[:1200]
processing_nums = int(sys.argv[1])-1
pool = multiprocessing.Pool()
for index, files in enumerate(data_iter(processing_nums, all_files)):
    if index == 14:
        print(f"exec process {index} len(files)={len(files)}")
        func(index, files,)
pool.close()
pool.join()
