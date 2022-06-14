# %%
import sys
import os
import pdb
from load_ncbi_taxinfo import *
from collections import Counter
from ete3 import NCBITaxa
ncbi = NCBITaxa()
# %%
print("load taxonomy done!")
f = open("/data/home/wlzhang/classfication/software/kraken2/sourmash/data_out/sourmash.OPAL", "w")
for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
    version = "1.0.0"
    sample_id = str(i)
    file_name = f"/data/home/wlzhang/classfication/software/kraken2/sourmash/data_out/{i}.fastq"
    print(f"handle {i}  ")
    f.write(f"@SampleID:{sample_id}\n")
    f.write(f"@Version:{version}\n")
    f.write(
        f"@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n")
    f.write(f"\n")
    f.write(f"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
    total = 0
    total_taxid = []
    for line in open(f"{file_name}.sourmash", "r"):
        line = line.strip().split("%")
        percentage = line[0]
        if ";" in line[1]:
            name = line[1].split("/data/home")[0].split(";")[-1].strip()
        else:
            name = line[1].split("/data/home")[0].strip().split(" ")[-1]
        taxid = ncbi.get_name_translator([name])[name][0]
        rank = tax_id_to_rank[taxid]
        taxpath = get_id_path(taxid, *taxonomy)
        taxpath_text = "|".join([str(i) for i in taxpath])
        taxpathsn_text = "|".join(
            [str(tax_id_to_name[i]) for i in taxpath])
        text = f"{taxid}\t{rank}\t{taxpath_text}\t{taxpathsn_text}\t{percentage}\n"
        f.write(text)

# %%
