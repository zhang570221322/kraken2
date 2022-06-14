
import sys,os,pdb
from load_ncbi_taxinfo import *
print("load taxonomy done!")
prefix="/data/home/wlzhang/classfication/software/kraken2/Mash/data_out"
for i in [0,1,2,3,4,5,6,7,8,9]:
# for i in [0]:
    print(f"handle {i} profilling")
    version = "1.0.0"
    sample_id=str(i)
    file_name = f"{prefix}/{i}.fastq.mash"
    with open(f"{file_name}.OPAL","w") as f:
        f.write(f"@SampleID:{sample_id}\n")
        f.write(f"@Version:{version}\n")
        f.write(f"@Ranks:superkingdom|phylum|class|order|family|genus|species|strain\n") 
        f.write(f"\n")
        f.write(f"@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for line in open(f"{file_name}","r"):
            if line[0]=="%":
                continue
            line=line.strip().split("\t")
            percentage=line[0]
            taxid = int(line[-2])
            rank = tax_id_to_rank[taxid]
            taxpath = get_id_path(taxid,*taxonomy)
            taxpath_text =  "|".join([str(i) for i in taxpath])
            taxpathsn_text = "|".join([str(tax_id_to_name[i]) for i in taxpath])
            text = f"{taxid}\t{rank}\t{taxpath_text}\t{taxpathsn_text}\t{percentage}\n"
            f.write(text)
