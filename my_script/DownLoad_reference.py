# %%
import multiprocessing
import os
# 进程数
processing_nums = 10

# %%
suffix = "_genomic.fna.gz"
if not os.path.exists("assembly_summary.txt"):
    print("Download bacteria/assembly_summary.txt")
    os.system(
        "wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")


lines = open("assembly_summary.txt", "r").readlines()[1:]


def data_iter(batch, data):
    """
    # 均分data
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]


def func(muti_lines):
    for line in muti_lines:
        temp = line.split("\t")
        if len(temp) > 3:
            url_prefix = temp[-4]
            file_name = url_prefix.split("/")[-1]+suffix
            download_url = "/".join([url_prefix,
                                     file_name])
            os.system("wget -q  -P reference " + download_url)
            os.system("gunzip  ./reference/"+file_name)


# 多进程
pool = multiprocessing.Pool()
for lines_pool in data_iter(processing_nums, lines):
    pool.apply_async(func, (lines_pool,))
pool.close()
pool.join()


# %%
