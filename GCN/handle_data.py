import pdb
import numpy as np
import scipy.sparse as sp
from torch.utils import data
import torch


def de_to_dense(adj):
    res = []
    for index in range(adj.shape[0]):
        res.append(adj[index].to_dense())
    return torch.stack(res, 0)


def csc_adj_to_torch_spares(test, train_size):
    res_train = []
    res_test = []
    for index, adj in enumerate(test):
        Acoo = adj.tocoo()
        Apt = torch.sparse.LongTensor(torch.LongTensor([Acoo.row.tolist(), Acoo.col.tolist()]),
                                      torch.LongTensor(Acoo.data.astype(np.int16)), torch.Size([1024, 1024]))
        if index < train_size:
            res_train.append(Apt)
        else:
            res_test.append(Apt)
    return torch.stack(res_train, 0), torch.stack(res_test, 0)


def de_csc_adj(test):
    res = []
    for adj in test:
        Apt = torch.from_numpy(adj.toarray())
        res.append(Apt)
    return torch.stack(res, 0)


def standardization(data):
    mu = np.mean(data, axis=0)
    sigma = np.std(data, axis=0)
    return (data - mu) / sigma


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


def load_data(preifx="/home/yegpu/zwl/data/8k_reads", data_size=8000):
    x = np.load(f"{preifx}/x.npy", allow_pickle=True)[:data_size]
    x = standardization(x)
    adjacent_matrixs = np.load(
        f"{preifx}/adjacent_matrixs.npy", allow_pickle=True)[:data_size]
    Y = np.load(f"{preifx}/y_genus_one_hot.npy", allow_pickle=True)[:data_size]
    return x, adjacent_matrixs, Y


def handle_data():
    X, adj, Y = load_data()
    X = torch.from_numpy(X)
    X = X.unsqueeze(-1)
    Y = torch.from_numpy(Y)
    return X, adj, Y


def load_array(data_arrays, batch_size, is_train=True):  # @save
    """构造一个PyTorch数据迭代器"""
    dataset = data.TensorDataset(*data_arrays)
    return data.DataLoader(dataset, batch_size, shuffle=is_train)


def getdata(train_ratio=0.8):
    X, adj, Y = handle_data()
    train_size = int(train_ratio*len(X))
    adj = de_csc_adj(adj)
    out_dim = Y.shape[-1]
    # Y = Y.max(axis=-1).indices
    # 用来训练
    train_data = (X[:train_size], adj[:train_size], Y[:train_size])
    # train_data = (X, adj, Y)
    # 用来测试
    test_data = (X[train_size:], adj[train_size:], Y[train_size:])
    return train_data, out_dim, test_data


def getdata2(train_ratio=0.8):
    X, adj, Y = handle_data()
    train_size = int(train_ratio*len(X))
    adj_train, adj_test = csc_adj_to_torch_spares(adj, train_size)
    # 用来训练
    train_data = (X[:train_size], adj_train, Y[:train_size])
    # 用来测试
    test_data = (X[train_size:], adj_test, Y[train_size:])
    return train_data, test_data
# out_dir = "/home/yegpu/zwl/data/2k_reads"
# x,adj,y = load_data()
# np.save(f"{out_dir}/x", x[:2000])
# np.save(f"{out_dir}/adjacent_matrixs", adj[:2000])
# np.save(f"{out_dir}/y_genus_one_hot", y[:2000])
