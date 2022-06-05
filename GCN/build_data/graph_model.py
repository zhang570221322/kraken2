import pdb
DEFAULT_K = 12

KMER_DIC = {
    "A": [1,0,0,0],
    "C": [0,1,0,0],
    "G": [0,0,1,0],
    "T": [0,0,0,1],
    "a": [1,0,0,0],
    "c": [0,1,0,0],
    "g": [0,0,1,0],
    "t": [0,0,0,1],

}
REVERSE = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "a": "T",
    "c": "G",
    "g": "C",
    "t": "C",

}
class Read:
    def __init__(self, _id, seq, annotation, quality, raw):
        self.id = _id
        self.seq = seq
        self.annotation = annotation
        self.quality = quality
        self.raw = raw
def feature_str(read):
    return [KMER_DIC[ch] for ch in read]

def reverse_sqe(read:Read):
    rev_content = read.seq
    read.seq =[REVERSE[i] for i in rev_content][::-1]
    # WARNNING quality
    return read

def fasta2fastq(read:Read,k=220):
    dp_length = len(read.seq)
    if(dp_length) < k:
        raise Exception(f"The sequence {read.id} dp < {k}!!")
    start = 0
    end = k
    while end <= dp_length:
        yield Read(read.id,read.seq[start:end],read.annotation,None,None)
        start += 1
        end += 1
def get_taxid_kraken2(seq_id):
    tax_id = seq_id.split("|kraken:taxid|")[-1].split("_")[0]
    if not tax_id.isdigit():
        tax_id = seq_id.split("|kraken:taxid|")[-1].split(" ")[0]
    return tax_id



class ReadGenerator:
    def __init__(self, file_name, seq_type="reads"):
        """[ReadGenerator]
        Args:
            file_name ([str]): [the file full name]
            seq_type (str): [reads(4rows),fasta]. Defaults to "reads".
        """
        self.file_name = file_name
        self.seq_type = seq_type

    def get_read_length(self):
        count = 0
        for _ in self.read_Generator():
            count += 1
        return count
    def get_fasta2read_length(self):
        count = 0
        times=0
        for read in self.read_Generator():
            times+=1
            count += len(read.seq)
        return count-times*119
    def read_Generator(self,default_fastq_length=220):
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
                    if read.seq != "":
                        read.id = raw[0][1:].strip()
                        read.seq = sequence_content
                        read.raw = "".join(raw[:-1])[:-1]
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
        elif seq_type =="fasta2fastq":
            raw = []
            sequence_content = ""
            read = Read(None, None, None, None, None)
            for line in open(self.file_name, "r"):
                raw.append(line)
                if line.startswith(">"):
                    if read.seq != None:
                        read.id = raw[0][1:].strip()
                        read.seq = sequence_content
                        read.raw = "".join(raw[:-1])[:-1]
                        for fastq_read in fasta2fastq(read,default_fastq_length):
                            yield fastq_read
                        for fastq_read in  fasta2fastq(reverse_sqe(read),default_fastq_length):
                            yield fastq_read
                        sequence_content = ""
                        raw = [line]
                        read = Read(None, None, None, None, None)
                    
                else:
                    sequence_content += line.strip()
            # 处理最后一个
            read.id = raw[0][1:].strip()
            read.seq = sequence_content
            read.raw = "".join(raw)
            for fastq_read in fasta2fastq(read,default_fastq_length):
                yield fastq_read
            for fastq_read in  fasta2fastq(reverse_sqe(read),default_fastq_length):
                yield fastq_read
    def Kmer_index_Generator(self, read, k=DEFAULT_K):
        dp_length = len(read.seq)
        seq_feature = feature_str(read.seq)
        if(dp_length) < k:
            raise Exception(f"The sequence {read.id} dp < {k}!!")
        start = 0
        end = k
        while end <= dp_length:
            yield read.seq[start:end], (start, end), seq_feature[start:end]
            start += 1
            end += 1


class FeatureSpace(dict):
    def __init__(self, k=DEFAULT_K) -> None:
        super().__init__()
        self.k = k

    def init(self):
        def generater(temp, k):
            if len(temp) == k:
                self["".join(temp)] = [0]
                return
            for base in "ATCG":
                temp.append(base)
                generater(temp, k)
                temp.pop()
        generater([], self.k)
        for index, kmer in enumerate(self.keys()):
            self[kmer] = [index, 0, 0]+feature_str(kmer)


feature_space_init = FeatureSpace()
if __name__ == "__main__":
    file_name = "./test.fna"
    read_generator = ReadGenerator(file_name, "fasta")
    for read in read_generator.read_Generator():
        print(read.id)
        print(read.raw)
