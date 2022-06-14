# %%
import multiprocessing
import os


# %%


def get_taxids(ancestor_taxid):
    ancestor_taxid = str(ancestor_taxid)
    parent_childrens = {}
    for line in open("/data/home/wlzhang/classfication/kraken2_db/2021_11_01/taxonomy/nodes.dmp", "r"):
        temp = line.split("\t|\t")
        childrens = parent_childrens.get(temp[1])
        if childrens:
            childrens.append(temp[0])
        else:
            parent_childrens[temp[1]] = [temp[0]]
    taxid = set()

    def dfs(childrens, parent):
        if childrens is None:
            taxid.add(parent)
            return
        for parent in childrens:
            taxid.add(parent)
            dfs(parent_childrens.get(parent), parent)
        # get taxid
        dfs(parent_childrens[ancestor_taxid], ancestor_taxid)
    return taxid


# 下载18601（细菌）等级下的所有物种序列
ancestor_taxid = 18601
taxid = get_taxids(ancestor_taxid)

suffix = "_genomic.fna.gz"
if not os.path.exists("assembly_summary.txt"):
    os.system(
        "wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")
total = 0
for line in open("assembly_summary.txt", "r"):
    temp = line.split("\t")
    if len(temp) > 3:
        if temp[5] in taxid:
            total += 1
print(f"ancestor_taxid:{ancestor_taxid},count:{total}")


lines = open("assembly_summary.txt", "r").readlines()[1:]


def data_iter(batch, data):
    """
    # 均分taxid
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]


def func(muti_taxid, muti_lines):
    for line in muti_lines:
        temp = line.split("\t")
        if len(temp) > 3:
            if temp[5] in muti_taxid:
                url_prefix = temp[-4]
                file_name = url_prefix.split("/")[-1]+suffix
                download_url = "/".join([url_prefix,
                                         file_name])
                os.system("wget -q  -P reference " + download_url)



# 多进程
thread= 10
pool = multiprocessing.Pool()
for muti_taxid in data_iter(thread, taxid):
    pool.apply_async(func, (muti_taxid, lines))
pool.close()
pool.join()

# %%
