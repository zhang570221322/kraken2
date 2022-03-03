import numpy as np
from ete3 import NCBITaxa
import sys
from sklearn import preprocessing
y_dir = sys.argv[1]


def get_genus(tax_id):
    loc = -1
    temp = ncbi.get_lineage(tax_id)
    target = temp[loc]
    while ncbi.get_rank([target])[target] != "genus":
        loc = loc-1
        target = temp[loc]
    return target


ncbi = NCBITaxa()
y = np.load(y_dir, allow_pickle=True)
for index, tax_id in enumerate(y):
    target = get_genus(tax_id)
    y[index] = target

y = y.reshape(-1, 1)
encoder = preprocessing.OneHotEncoder()
encoder.fit(y)
y = encoder.transform(y).toarray()
out_dir = f"{y_dir.split('.')[0]}_genus_one_hot"
np.save(out_dir, y)
