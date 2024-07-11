#!/usr/bin/env python


import argparse
import sys

from datetime import datetime
from pathlib import Path

from src.kmc import (count_kmers, create_input_file, create_kmer_histogram, 
                     calculate_cutoffs, dump_kmer_counts, get_hetkmers)
from src.utils import check_run


def parse_arguments():
    desc = "Calculate kmers occurrences from fastq/fasta files"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) TSV file containing accession, kind 
                    and files separated by commas, for example
                    ACC1    file1,file2
                    ACC2    file3
                    If there is more than one file for accession
                    kmers while be combined by union-sum'''
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output dir"
    parser.add_argument("--output_dir", "-o", type=str,
                        help=help_output, required=True)
    
    help_kmer_length = "(Required) kmer length"
    parser.add_argument("--kmer_length", "-k", type=int,
                        help=help_kmer_length, required=True)
    
    help_keep_intermediate = "(Optional) remove intermediate meryl files"
    parser.add_argument("--keep_intermediate", "-s", action="store_true",
                        default=False, help=help_keep_intermediate)
  
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = {}
    input_sequence = Path(parser.input_file)
    with open(input_sequence) as input_fhand:
        for line in input_fhand:
            if line:
                line = line.rstrip().split()
                name = line[0]
                files = line[1].split(",")
                if name in inputs:
                    raise ValueError("{} is duplicated".format(name))
                else:
                    inputs[name] = files        
    return {"inputs": inputs,
            "output": Path(parser.output_dir),
            "kmer_length": parser.kmer_length,
            "keep_intermediate": parser.keep_intermediate}


def main():
    arguments = get_arguments()
    logdate = "Kmer_counting_"+datetime.now().strftime("%d_%m_%Y-%H_%M_%S") + ".log"
    log_fname = arguments["output"] / logdate
    if not arguments["output"].exists():
        arguments["output"].mkdir(exist_ok=True, parents=True)
        
    with open(log_fname, "w") as log_fhand:
        log_fhand.write("#Command Used: "+ "".join(sys.argv)+"\n")
        for name, input_files in arguments["inputs"].items():
            input_file_path = create_input_file(input_files, name, 
                                                arguments["output"])
            results = count_kmers(input_file_path, name, 
                                  arguments["output"], kmer_size=21, 
                                  threads=40, max_ram=64)
            log_fhand.write(check_run(results)+"\n")
            results = create_kmer_histogram(results["out_fpath"], name)
            log_fhand.write(check_run(results)+"\n")
            lower_bound, upper_bound = calculate_cutoffs(results["out_fpath"])
            results = dump_kmer_counts(Path(arguments["output"]), name, threads=40, 
                                       lower_bound=lower_bound, upper_bound=upper_bound)
            log_fhand.write(check_run(results)+"\n")
            results = get_hetkmers(Path(results["out_fpath"]), name)
            log_fhand.write(check_run(results)+"\n")
            exit()
    # counts = count_kmers(arguments)
    # with open(log_fname, "w") as log_fhand:
    #     )
    #     for name, results in counts.items():
    #         check = check_run(results)
    #         log_fhand.write("{}\t{}\n".format(name, check))



if __name__ == "__main__":
    main()