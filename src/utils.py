from subprocess import run


MAGIC_NUMS_COMPRESSED = [b'\x1f\x8b\x08', b'\x42\x5a\x68', 
                         b'\x50\x4b\x03\x04']

class UnionFind:
    def __init__(self, elementos):
        self.padre = {elemento: elemento for elemento in elementos}
        self.rango = {elemento: 0 for elemento in elementos}

    def find(self, elemento):
        if self.padre[elemento] != elemento:
            self.padre[elemento] = self.find(self.padre[elemento])
        return self.padre[elemento]

    def join(self, elemento1, elemento2):
        raiz1 = self.find(elemento1)
        raiz2 = self.find(elemento2)
        
        if raiz1 != raiz2:
            if self.rango[raiz1] > self.rango[raiz2]:
                self.padre[raiz2] = raiz1
            elif self.rango[raiz1] < self.rango[raiz2]:
                self.padre[raiz1] = raiz2
            else:
                self.padre[raiz2] = raiz1
                self.rango[raiz1] += 1

def get_universe_size(fpaths):
    if len(fpaths) == 1:
        cmd = "wc -l {}".format(" ".join(fpaths))
    else:
        cmd = "cut -f1 {}| sort |uniq|wc -l".format(" ".join(fpaths))
    results = run(cmd, capture_output=True, shell=True)
    return int(results.stdout.decode().strip().split()[0])


def check_run(results):
    if results["returncode"] == 0:
        return "#SUCCESS: {}".format(results["command"])
    elif results["returncode"] == 99:
        return "#ALREADY_DONE: {}".format(results["command"])
    else:
        return "#FAIL: {} {}".format(results["command"], results["msg"])
    

def file_is_compressed(path):
    max_len = max(len(x) for x in MAGIC_NUMS_COMPRESSED)
    with open(path, 'rb') as fhand:
        file_start = fhand.read(max_len)
        fhand.close()
    for magic_num in MAGIC_NUMS_COMPRESSED:
        if file_start.startswith(magic_num):
            return True
    return False


def sequence_kind(path):
    #Checks if sequence file is in fasta or fastq format
    if file_is_compressed(path):
       cat_command = ["zcat"]
    else:
        cat_command = ["cat"]
    cat_command += [str(path), "|", "head", "-n 4"]
    run_cat_command = run("\t".join(cat_command), 
                          shell=True, capture_output=True)
    
    #fastq files headers starts with "@"
    #fasta file headers starts with (">")
    if run_cat_command.stdout.startswith(b'@'):
        return "fastq"
    elif run_cat_command.stdout.startswith(b'>'):
        return "fasta"
    else:
        msg = "{} doesn't seems to be a valid"
        msg += " fasta/fastq file"
        raise RuntimeError(msg.format(path.name))
    

def merge_dump_files(filepaths, out_filepath):
    filepaths = [str(filepath) for filepath in filepaths]
    out_filepath = "{}.dump".format(out_filepath)
    cmd = ["cat"] + filepaths 
    run_ = run(" ".join(cmd), shell=True, capture_output=True)
    counts = {}
    for line in run_.stdout.decode().split("\n"):
        if line:
            kmer, count = line.rstrip().split()
            if kmer not in counts:
                counts[kmer] = int(count)
            else:
                counts[kmer] += int(count)
    with open(out_filepath, "w") as out_fhand:
        for kmer, count in counts.items():
            out_fhand.write("{}\t{}\n".format(kmer, count)) 
    return out_filepath


def merge_hetkmer_files(filepaths, out_fdir):
    out_fpath = out_fdir / "hetkmers_sequences_merged.tsv"
    cmd = "sort -u {} > {}".format(" ".join(filepaths), str(out_fpath))
    print(cmd)
    run(cmd, shell=True)

def get_kmer_value(filepath, kmer):
    cmd = "grep -w \"{}\" {} | head -1".format(kmer, str(filepath))
    value = run(cmd, shell=True, capture_output=True)
    if not value.stdout:
        return 0
    else:
        return int(value.stdout.decode().rstrip().split()[1])
    
def merge_temporary_files(tmp_dir, out_dir, suffix):
    if suffix == ".estim":
        out_fpath = out_dir / "kmer_estimators.tsv"
    if suffix == ".index":
        out_fpath = out_dir / "kmer_index.tsv"
    cmd = "ls -v {}/*{}|xargs cat > {}".format(str(tmp_dir), suffix, str(out_fpath))
    run(cmd, shell=True)
