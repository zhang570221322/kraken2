import os
import uuid
import sys

file_dir = sys.argv[1]
out_dir = str(uuid.uuid1())
with open(out_dir, "w") as f:
    for line in open(file_dir, "r"):
        if line.startswith("@") or line.startswith(">"):
            f.write(line)
        else:
            line = line.upper()
            if "N" in line:
                line = line.replace("N", "")
            if line.strip() == "":
                continue
            f.write(line)
os.system(f"rm {file_dir}")
os.system(f"mv {out_dir} {file_dir}")
