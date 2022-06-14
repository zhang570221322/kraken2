import pdb
import os
import sys
import random
reference_data = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/Bacteria_artificial/GCF_000019.fna"
wgsim = "/data/home/wlzhang/classfication/software/stimulation/wgsim/wgsim"
out_dir = sys.argv[1]+"/"


def create(file):
    error = round(random.uniform(0.02, 0.18), 2)
    os.system(
        f"{wgsim} -N50000 -1220 -e {error}  {reference_data} {file} /dev/null >/dev/null")
    # 生成金标准文件
    with open(f"{file}.binning", "w") as f:
        for line in open(f'{file}', "r"):
            if line[0] == ("@"):
                read_id = line.strip()
                read_taxid = read_id.split("|")[-1].split("_")[0]
                f.write(f"{read_id[1:]}\t{read_taxid}\n")


files = [str(i)+".fastq" for i in range(10)]
for i in files:
    create(out_dir+i)
