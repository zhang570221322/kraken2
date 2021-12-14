import numpy as np
import scipy.sparse as sp
from torch.utils import data

import torch


def normalize_feature(features):
    """
    Row-normalize feature matrix and convert to tuple representation
    """
    rowsum = np.array(features.sum(
        2), dtype=np.float64)  # get sum of each row, [2708, 1]
    r_inv = np.power(rowsum, -1)   # 1/rowsum, [2708]
    r_inv[np.isinf(r_inv)] = 0.  # zero inf data
    for i in range(r_inv.shape[1]):
        r_mat_inv = sp.diags(r_inv[:, i])
        features[:, i, :] = r_mat_inv.dot(features[:, i, :])
    return features


def load_data(N=220):
    preifx = "/data/home/wlzhang/classfication/kraken2_db/2021_11_01/library/stimulation/test/dp/data"
    x1 = list(np.load(f"{preifx}/x1.npy", allow_pickle=True))
    x2 = np.asarray(np.load(f"{preifx}/x2.npy",
                    allow_pickle=True), dtype=np.float32)
    adjacent_matrixs = np.load(
        f"{preifx}/adjacent_matrixs.npy", allow_pickle=True)
    Y = list(np.load(f"{preifx}/Y.npy", allow_pickle=True))
    # 还原x1,x2至x
    for index, temp in enumerate(zip(x1, x2)):
        temp_x1, temp_x2 = temp
        need = list(np.zeros((N-len(temp_x1),), dtype=np.float32))
        x1[index] = [x1[index]+need, temp_x2]
    X = np.asarray(x1, dtype=np.float32)
    # 还原邻接矩阵
    adjacent_matrixs_np = np.zeros((X.shape[0], N, N), dtype=np.int8)
    for index, ab in enumerate(adjacent_matrixs):
        for a, b in ab:
            adjacent_matrixs_np[index, a, b] = 1
    # 还原标签Y
    for index, temp in enumerate(Y):
        need = list(np.zeros((N-len(temp),), dtype=np.float32))
        Y[index] = Y[index]+need
    Y = np.asarray(Y, dtype=np.float32)
    return X, adjacent_matrixs_np, Y


def handle_data():
    X, adjacent_matrixs, Y = load_data()
    X = normalize_feature(X)
    X = X.transpose(0, 2, 1)

    X = torch.from_numpy(X)
    adj = torch.from_numpy(adjacent_matrixs)
    Y = torch.from_numpy(Y)
    Y = Y.unsqueeze(-1)
    return X, adj, Y


def load_array(data_arrays, batch_size, is_train=True):  # @save
    """构造一个PyTorch数据迭代器"""
    dataset = data.TensorDataset(data_arrays)
    return data.DataLoader(dataset, batch_size, shuffle=is_train)


X, adj, Y = handle_data()
train_size = int(0.8*len(X))
# 用来训练
train_data = (X[:train_size], adj[:train_size], Y[:train_size])
# 用来测试
test_data = (X[train_size:], adj[train_size:], Y[train_size:])
# train_data_iter = load_array(train_data,  train_size//100)


# x = torch.from_numpy(X[:10, ].copy())
# supports = torch.from_numpy(adjacent_matrixs[:10, ].copy())
# y = torch.from_numpy(Y[:10, ].copy())
# torch.save([x, supports, y], 'test-data')
# device = torch.device('cpu')
