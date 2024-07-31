import pandas as pd

from io import StringIO
from subprocess import run
from pathlib import Path
from shutil import rmtree as remove_dir

from src.utils import reformat_lines


###Operation for getting kmers from paired end reads
'''meryl union-sum [count k=21 SRA_R1.fastq.gz output SRA_R1.count] 
    [count k=21 SRA_R2.fastq.gz output SRA_R2.count]  
    output SRA.union_sum'''


def _run_meryl(name, fpaths, kmers_outfpath, kmer_length, keep_intermediate=False):
    out_fpath = "{}.count".format(str(kmers_outfpath/name))     


    if len(fpaths) == 1:
         fpath = "".join(fpaths)
         cmd = "meryl count k={} {} output {}".format(kmer_length, 
                                                      "".join(fpath),
                                                      out_fpath)
    else:
        bin = ["meryl union-sum"]
        for fpath in fpaths:
            count = ["[count k={} {} output {}.countpart]".format(kmer_length, 
                                                                  fpath,
                                                                  kmers_outfpath/Path(fpath).stem)]
            bin += count

        bin += ["output", out_fpath]

        cmd = " ".join(bin)

    if Path(out_fpath).exists():
         results = {"command": cmd, "returncode": 99,
                    "msg": "Output file already exists", "out_fpath": out_fpath}
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        results = {"command": cmd, "returncode": run_.returncode, "name": name,
                    "msg": run_.stderr.decode(), "out_fpath": out_fpath}
    if not keep_intermediate:
         for intermediate_file in kmers_outfpath.glob("*.countpart"):
              remove_dir(intermediate_file)
    
    return results


def count_kmers(arguments):
    results = {}
    for name, files in arguments["inputs"].items():
            out = _run_meryl(name, files, arguments["output"], arguments["kmer_length"], 
                             keep_intermediate=arguments["keep_intermediate"])
            results[name] = out
    return results


def get_kmer_counts_dataframe(filepath, hetkmers={}):
    name = str(filepath.stem).split(".")[0]
    count_colname = "{}_kmer_count".format(name)
    with open(filepath) as filehand:
         lines = [[line.split()[0], int(line.rstrip().split()[1])] for line in filehand]
    lines = reformat_lines(lines, hetkmers=hetkmers)
    io = StringIO("\n".join(lines)) 
    df = pd.read_csv(io, delimiter="\t", names=["kmer", "COUNT"])
    df = df.groupby(["kmer"]).COUNT.sum().reset_index()
    df.columns = ["kmer", count_colname]
    df = df.sort_values(by=["{}_kmer_count".format(name)], ascending=False)
    return df