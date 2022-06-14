import os
import sys
import multiprocessing


def data_iter(batch, data):
    """
    # 均分data
    """
    data = list(data)
    num_examples = len(data)
    step_size = num_examples//batch
    for i in range(0, num_examples, step_size):
        yield data[i: min(i + step_size, num_examples)]


def func(index, files):
    for file_name in files:
        file_dir = f"{prefix}/{file_name}"
        out_file = f"{our_dir}/{file_name}"
        with open(out_file, "w") as f:
            for line in open(file_dir, "r"):
                if line.startswith(">"):
                    f.write(line)
                else:
                    line = line.upper()
                    if "N" in line:
                        line = line.replace("N", "")
                    f.write(line)


prefix = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/Bacteria_artificial/reference"
our_dir = "/data/home/wlzhang/classfication/software/kraken2/test/stimulation_test/Bacteria_artificial/reference_delete_n"
all_files = os.listdir(prefix)
processing_nums = int(sys.argv[1])
pool = multiprocessing.Pool()
for index, files in enumerate(data_iter(processing_nums, all_files)):
    print(f"exec process {index} len(files)={len(files)}")
    pool.apply_async(func, (index, files,))
pool.close()
pool.join()
