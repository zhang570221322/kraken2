import numpy as np
import torch
import util as util


def load_data(N=220):
    preifx = "/data/home/wlzhang/classfication/kraken2_db/2021_11_01/library/stimulation/test/dp/data"
    x1 = list(np.load(f"{preifx}/x1.npy", allow_pickle=True))
    x2 = list(np.load(f"{preifx}/x2.npy", allow_pickle=True))
    adjacent_matrixs = np.load(
        f"{preifx}/adjacent_matrixs.npy", allow_pickle=True)
    Y = list(np.load(f"{preifx}/Y.npy", allow_pickle=True))
    # 还原x1,x2至x
    for index, temp in enumerate(zip(x1, x2)):
        temp_x1, temp_x2 = temp
        need = list(np.zeros((N-len(temp_x1),), dtype=np.float64))
        x1[index] = [x1[index]+need, temp_x2]
    x = np.asarray(x1, dtype=np.float64)
    # 还原邻接矩阵
    adjacent_matrixs_np = np.zeros((x.shape[0], N, N), dtype=np.int8)
    for index, ab in enumerate(adjacent_matrixs):
        for a, b in ab:
            adjacent_matrixs_np[index, a, b] = 1
    # 还原标签Y
    for index, temp in enumerate(Y):
        need = list(np.zeros((N-len(temp),), dtype=np.float64))
        Y[index] = Y[index]+need
    Y = np.asarray(Y, dtype=np.float64)
    return x, adjacent_matrixs, Y


X, adjacent_matrixs, Y = load_data()

features = util.preprocess_features(X)  # [49216, 2], [49216], [2708, 1433]
supports = util.preprocess_adj(adjacent_matrixs)

device = torch.device('cpu')
i = torch.from_numpy(features[0]).long().to(device)
v = torch.from_numpy(features[1]).to(device)
feature = torch.sparse.FloatTensor(
    i.t(), v, features[2]).float().to(device)

i = torch.from_numpy(supports[0]).long().to(device)
v = torch.from_numpy(supports[1]).to(device)
support = torch.sparse.FloatTensor(
    i.t(), v, supports[2]).float().to(device)
