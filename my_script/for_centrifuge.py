import os


prefix = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/Bacteria_artificial/reference"
our_dir = "/data/home/wlzhang/classfication/software/kraken2/centrifuge/db"
all_files = os.listdir(prefix)
out_file = f"{our_dir}/Bacteria_artificial.fna"
fna = open(out_file, "w")
nc_taxid = open(f"{out_file}.binning", "w")
for file_name in all_files:
    file_dir = f"{prefix}/{file_name}"
    for line in open(file_dir, "r"):
        if line.startswith(">"):
            line = line.split("|kraken:taxid|")
            taxid = line[1].split(" ")[0]
            new_line = f"{line[0]}\n"
            nc_taxid.write(f"{line[0][1:]}\t{taxid}\n")
            fna.write(new_line)
        else:
            fna.write(line)
