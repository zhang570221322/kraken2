#%%
from load_ncbi_taxinfo import *
taxid_path,name_path,rank_path=get_id_path(1496,*taxonomy)
print("cat#### ")
print(get_id_path(9685, *taxonomy))
print("dog#### ")
print(get_id_path(9615, *taxonomy))