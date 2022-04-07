from autogl.datasets import utils
from ReadsDataSet_kraken_all import ReadsDataSet
import torch, numpy as np ,random
from autogl.solver import AutoGraphClassifier 
import warnings
import pdb
warnings.filterwarnings("ignore")
print("start to load data...")
dataset = ReadsDataSet(root='./')
print("load data done!")

seed=12345
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
utils.graph_random_splits(dataset, train_ratio=0.8, val_ratio=0.1, seed=seed)
solver = AutoGraphClassifier ()
solver.fit(
    dataset,
    time_limit=3600,
    train_split=0.8,
    val_split=0.05,
)
solver.get_leaderboard().show()
print('best single model:\n', solver.get_leaderboard().get_best_model(0))
from autogl.module.train import Acc
predicted = solver.predict_proba()
print('Test accuracy: ',Acc.evaluate(predicted, dataset.data.y[dataset.test_index].cpu().detach().numpy()))
