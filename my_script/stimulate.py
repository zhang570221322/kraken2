import pyfastx
import os
import random
import uuid
prefix = "/data/home/wlzhang/classfication/"
dir_path = f"{prefix}kraken2_db/2021_11_01/taxonomy/Bacteria_artificial/reference"
a = os.listdir(dir_path)
kraken2_db = f"{prefix}kraken2_db/2021_11_01/"
u = "1"
for file in random.sample(a, 2000):
    complete_file = f"{dir_path}/{file}"
    os.system(
        f"/data/home/wlzhang/classfication/software/stimulation/wgsim/wgsim -N9000 -1250  {complete_file} temp_{file}.fastq /dev/null > /dev/null")
    os.system(f"cat temp_{file}.fastq  >> {u}.fastq")
    os.remove(f"temp_{file}.fastq")
# 生成金标准文件
with open(f"{u}.fastq.binning", "w") as f:
    for read in pyfastx.Fastq(f'{u}.fastq'):
        read_id = read.name
        read_taxid = read.name.split("|")[-1].split("_")[0]
        f.write(f"{read_id}\t{read_taxid}\n")
