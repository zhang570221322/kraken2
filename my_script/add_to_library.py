# %%
import os
prefix = "/data/home/wlzhang/classfication/"
dir_path = f"{prefix}kraken2_db/2021_11_01/taxonomy/Bacteria_artificial/reference"
a = os.listdir(dir_path)
kraken2_build = f"{prefix}software/kraken2/kraken2_weight/kraken2_build/kraken2-build"
kraken2_db = f"{prefix}kraken2_db/2021_11_01/"
# %%
for file in a:
    file = f"{dir_path}/{file}"
    os.system(
        f"{kraken2_build} -t 10 --no-masking --add-to-library {file} --db {kraken2_db}")
