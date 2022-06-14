import sys
base_file = sys.argv[1]

f = open(f"{base_file}.diff.output", "w")
# f2 = open(f"{base_file}.same.output", "w")
_base = {}
for line in open(f"{base_file}.origin.output", "r"):
    temp = line.split("\t")
    _base[temp[1]] = [temp[2], line]
for line in open(f"{base_file}.weight.output", "r"):
    temp = line.split("\t")
    if _base[temp[1]][0] != temp[2]:
        # f.write(_base[temp[1]][1])
        f.write(line)
    # else:
    #     f2.write(line)
f.close()
# f2.close()
