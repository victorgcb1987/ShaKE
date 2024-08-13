#!/usr/bin/env python

import argparse
import sys

from pathlib import Path
from shutil import rmtree as remove_dir
from src.kmer import calculate_kmer_estimators
from src.utils import merge_temporary_files


def parse_arguments():
    desc = "Create a table with "
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) Directory with meryl counts.
                    Subolders have the suffix .count'''
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output dir"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output, required=True)
    
    help_output = "(Required) Dump file with universe of kmers"
    parser.add_argument("--universe", "-u", required=True,
                        type=str)

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = sorted([dir for dir in Path(parser.input_dir).glob("*.dump")])
    return {"inputs": inputs,
            "output": Path(parser.output),
            "universe": Path(parser.universe),
            "N": len(inputs)}


def main():
    arguments = get_arguments()
    tmp_dir = arguments["output"] / "tmp"
    tmp_dir.mkdir(parents=False, exist_ok=True)
    with open(tmp_dir / "Header.estim", "w") as header_fhand:
        header_fhand.write("Kmer\tDiversity\tSpecifity\n")
    
    with open(arguments["universe"]) as universe_fhand:
        k = 0
        for line in universe_fhand:
            k += 1
            kmer = line.split()[0]
            with open(tmp_dir / "K{}.index".format(k), "w") as index_fhand:
                index_fhand.write("K{}\t{}\n".format(k, kmer))
            with open(tmp_dir / "K{}.estim".format(k), "w") as estimator_fhand:
                diversity, specifity = calculate_kmer_estimators(arguments["inputs"], arguments["N"], kmer)
                estimator_fhand.write("K{}\t{}\t{}\n".format(k, diversity, specifity))
        merge_temporary_files(tmp_dir, arguments["output"], ".estim")
        merge_temporary_files(tmp_dir, arguments["output"], ".index")
        remove_dir(tmp_dir)
        

if __name__ == "__main__":
    main()