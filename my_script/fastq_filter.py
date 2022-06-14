import pdb
import sys
from pojo import ReadGenerator
filter1 = [line.strip() for line in open(sys.argv[2], "r")]
file = sys.argv[1]


right = open(f"{file}.right", "w")
error = open(f"{file}.error", "w")
for read in ReadGenerator(f"{file}").get_Generator():
    if read.id in filter1:
        error.write(read.raw)
    else:
        right.write(read.raw)
