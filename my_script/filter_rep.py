import sys
import os
from pojo import ReadGenerator

fastq = sys.argv[1]
fastq_binning = fastq+".binning"
fastq_new = fastq+".new"
fastq_binning_new = fastq_binning+".new"
_fastq_id = set()

f1 = open(fastq_new, "w")
f2 = open(fastq_binning_new, "w")
for read in ReadGenerator(fastq).get_Generator():
    if read.id in _fastq_id:
        continue
    f1.write(read.raw)
    _fastq_id.add(read.id)
for line in open(fastq_binning, "r"):
    seq_id, taxo = line.strip().split("\t")
    if seq_id in _fastq_id:
        f2.write(line)
os.system(f"rm {fastq}")
os.system(f"rm {fastq_binning}")
os.system(f"mv {fastq_new} {fastq}")
os.system(f"mv {fastq_binning_new} {fastq_binning}")
