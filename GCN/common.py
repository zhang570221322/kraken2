# %%
import torch
import matplotlib.pyplot as plt
import arrow
from time import time
import functools
from torch_geometric.loader import DataLoader
from tqdm import tqdm
import pdb

def get_tax_list(taxs, target,mode=1):
    if mode==1:
        pre = taxs["root"]
        res = []
        try:
            for i in target:
                i = int(i)
                if i == 0:
                    return res
                tax = taxs[pre][i-1]
                res.append(tax)
                pre = tax
        except KeyError as err:
            return res
        except IndexError as err :
            return res
        return res
    elif mode == 2:
        res = []
        for first,second in enumerate(target):
            key = f"{first}_{second}"
            res.append(taxs[key])
        return res

def AUTO_INCREMENT(n=0):
    while True:
        yield n
        n = n+1


def standardization(data):
    mu = torch.mean(data, axis=0)
    sigma = torch.std(data, axis=0)
    return (data - mu) / sigma


def DataLoaderTqdm(data, *arg, **kwarg):
    length = len(data)
    batch_size = kwarg["batch_size"]
    # res = tqdm(DataLoader(data, *arg, **kwarg),
    #            total=length//batch_size, unit=f"batch({batch_size})")
    res = DataLoader(data, *arg, **kwarg)
    return res


def date_print(text):
    print(f"{date_str()}:{text}")


def date_str(format="YYYY_MM_DD_HH_mm_ss"):
    return arrow.now().format(format)


def log_time(text=""):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kw):
            # 定义开始时间
            start = time()
            res = func(*args, **kw)
            # 定义结束时间
            end = time()
            print(f'{text} {func.__name__}() execute time: {end-start:.2f}s')
            return res
        return wrapper
    return decorator


class Arg:
    def __init__(self, num_epochs, learning_rate, weight_decay, batch_size) -> None:
        self.num_epochs = num_epochs
        self.lr = self.learning_rate = learning_rate
        self.weight_decay = weight_decay
        self.batch_size = batch_size


def my_plot(*args, **kwargs, ):

    arg = kwargs.get("arg")
    args = list(args)
    titles = args.pop()
    nums_plot = len(args)
    figsize = (14, 4*nums_plot)
    plt.figure(figsize=figsize)
    x = list(range(len(args[0])))
    if arg:
        text = f"num_epochs:{arg.num_epochs}\nlearning_rate:{arg.learning_rate}\nweight_decay:{arg.weight_decay}\nbatch_size:{arg.batch_size}"
        plt.figtext(.00, .01,  text, fontsize=15)
    temp = nums_plot*100+11
    for y, title in zip(args, titles):
        plt.subplot(temp)
        plt.plot(x, y, color='r', marker='.', linestyle='-')
        plt.title(title)
        temp += 1
    file_name = date_str()
    plt.savefig(f"{file_name}.jpg")


# %%
if __name__ == "__main__":
    # test
    loss_record = train_acc_record = val_acc_record = test_acc_record = [
        1, 2, 3, 4, 5]
    arg = Arg(20, 0.01, 0.75, 1000)
    my_plot(loss_record, train_acc_record,
            val_acc_record, test_acc_record, ["Loss", "Train Auc", "Val Auc", "Test Auc"], arg=arg)

# %%
