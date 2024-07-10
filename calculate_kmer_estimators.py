#!/usr/bin/env python

import argparse
import sys

from datetime import datetime
from pathlib import Path

from src.kmer import get_kmer_counts_dataframe
from src.dataframes import (merge_kmer_dataframes, calculate_kmer_count_diversity, 
                            calculate_kmer_count_specifity, calculate_sample_kmer_diversity,
                            calculate_sample_kmer_specifity)


def parse_arguments():
    desc = "Create a table with "
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = '''(Required) Directory with meryl counts.
                    Subolders have the suffix .count'''
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output file"
    parser.add_argument("--output", "-o", type=str,
                        help=help_output, required=True)
    

    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    inputs = [dir for dir in Path(parser.input_dir).glob("*.count")]
    return {"inputs": inputs,
            "output": parser.output}


def main():
    arguments = get_arguments()
    dataframes = [get_kmer_counts_dataframe(dir) for dir in arguments["inputs"]]
    merged_dataframe_by_kmer = merge_kmer_dataframes(dataframes)
    print(merged_dataframe_by_kmer)
    merged_dataframe_by_sample = merged_dataframe_by_kmer.T.iloc
    header =  merged_dataframe_by_sample[0]
    merged_dataframe_by_sample = merged_dataframe_by_sample[1:]
    merged_dataframe_by_sample.columns = header
    calculate_sample_kmer_diversity(merged_dataframe_by_sample)
    calculate_sample_kmer_specifity(merged_dataframe_by_sample)
    merged_dataframe_by_sample.to_csv("test.csv", sep="\t")
    #calculate_sample_kmer_specifity(merged_dataframe_by_sample)
    # calculate_kmer_count_diversity(merged_dataframe_by_kmer)
    # calculate_kmer_count_specifity(merged_dataframe_by_kmer)
    # print(merged_dataframe_by_kmer)

if __name__ == "__main__":
    main()