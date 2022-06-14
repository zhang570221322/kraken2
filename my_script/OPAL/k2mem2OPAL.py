
import sys
import os
import pdb
from load_ncbi_taxinfo import *
from collections import Counter
print("load taxonomy done!")
f = open(f"/data/home/wlzhang/classfication/software/kraken2/K2Mem/data_out/k2mem.OPAL", "w")
for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
    version = "1.0.0"
    sample_id = str(i)
    file_name = f"/data/home/wlzhang/classfication/software/kraken2/K2Mem/data_out/{i}.fastq"
    type = "k2mem"
    print(f"handle {i} {type}")
    f.write(f"@SampleID:{sample_id}\n")
    f.write(f"@Version:{version}\n")
    f.write(
        f"@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
    f.write(f"\n")
    f.write(f"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
    total = 0
    total_taxid = []
    for line in open(f"{file_name}.{type}.output", "r"):
        line = line.strip().split("\t")
        taxid = int(line[2])
        total += 1
        taxid_path = get_id_path(taxid, *taxonomy)
        total_taxid += taxid_path
    temp = dict(Counter(total_taxid))
    temp.pop('')
    for taxid, times in temp.items():
        rank = tax_id_to_rank[taxid]
        taxpath = get_id_path(taxid, *taxonomy)
        taxpath_text = "|".join([str(i) for i in taxpath])
        taxpathsn_text = "|".join(
            [str(tax_id_to_name[i]) for i in taxpath])
        percentage = str(round(times/total, 10))
        text = f"{taxid}\t{rank}\t{taxpath_text}\t{taxpathsn_text}\t{percentage}\n"
        f.write(text)
