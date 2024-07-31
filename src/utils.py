from subprocess import run


MAGIC_NUMS_COMPRESSED = [b'\x1f\x8b\x08', b'\x42\x5a\x68', 
                         b'\x50\x4b\x03\x04']

def reformat_lines(kmer_count, hetkmers={}):
    return ["{}\t{}".format(hetkmers.get(count[0], count[0]), count[1]) for count in kmer_count]


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