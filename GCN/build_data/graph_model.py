
import pdb
import pickle
DEFAULT_K = 5


class Read:
    def __init__(self, _id, seq, annotation, quality, raw):
        self.id = _id
        self.seq = seq
        self.annotation = annotation
        self.quality = quality
        self.raw = raw


class ReadGenerator:
    def __init__(self, file_name, seq_type="reads") -> None:
        """[ReadGenerator]
        Args:
            file_name ([str]): [the file full name]
            seq_type (str): [reads(4rows),fasta]. Defaults to "reads".
        """
        self.file_name = file_name
        self.seq_type = seq_type

    def read_Generator(self):
        seq_type = self.seq_type
        if seq_type == "reads":
            step = 4
            cur = 0
            raw = []
            for line in open(self.file_name, "r"):
                cur += 1
                raw.append(line)
                if cur == step:
                    read = Read(raw[0][1:].strip(), raw[1].strip(),
                                raw[2].strip(), raw[3].strip(), "".join(raw))
                    cur = 0
                    raw = []
                    yield read
        elif seq_type == "fasta":
            raw = []
            sequence_content = ""
            read = Read(None, None, None, None, None)
            for line in open(self.file_name, "r"):
                raw.append(line)
                if line.startswith(">"):

                    read.id = line[1:].strip()
                    read.seq = sequence_content
                    read.raw = "".join(raw[:-1])[:-1]
                    if read.seq != "":
                        yield read
                        sequence_content = ""
                        raw = [line]
                        read = Read(None, None, None, None, None)
                else:
                    sequence_content += line.strip()
            # 处理最后一个
            read.id = raw[0][1:].strip()
            read.seq = sequence_content
            read.raw = "".join(raw)
            yield read

    def Kmer_index_Generator(self, read, k=DEFAULT_K):
        dp_length = len(read.seq)
        if(dp_length) < 8:
            raise Exception(f"The sequence {read.id} dp < {k}!!")
        start = 0
        end = k
        while end <= dp_length:
            yield read.seq[start:end], (start, end)
            start += 1
            end += 1


class FeatureSpace(dict):
    def __init__(self, k=DEFAULT_K) -> None:
        super().__init__()
        self.k = k

        def generater(temp, k):
            if len(temp) == k:
                self["".join(temp)] = [0]
                return
            for base in "ATCG":
                temp.append(base)
                generater(temp, k)
                temp.pop()
        generater([], self.k)
        for index, k in enumerate(self.keys()):
            self[k] = [index, 0]


feature_space_init = FeatureSpace()
# test
# file_name = "./test.fna"
# read_generator = ReadGenerator(file_name, "fasta")
# for read in read_generator.read_Generator():
#     print(read.id)
#     print(read.raw)
