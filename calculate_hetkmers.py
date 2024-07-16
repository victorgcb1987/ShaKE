#!/usr/bin/env python


import argparse
import sys

from datetime import datetime
from pathlib import Path

from src.kmc import calculate_hetkmers
from src.utils import concat_dump_files, check_run


def parse_arguments():
    desc = "get hetkmers between one or more kmc count dumps"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_input = "(Required) directory with dumps. dumps should have .dump filename extension"
    parser.add_argument("--input_file", "-i", type=str,
                        help=help_input, required=True)
    
    help_output = "(Required) output file"
    parser.add_argument("--output_dir", "-o", type=str,
                        help=help_output, required=True)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    input_fdir = Path(parser.input_file)
    output_fpath = Path(parser.output_dir)

    return {"input_folder": input_fdir,
            "out_fpath": output_fpath}


def main():
    arguments = get_arguments()
    files = [dump_file for dump_file in arguments["input_folder"].glob("*.dump")]
    combined_dump = concat_dump_files(files, arguments["out_fpath"])
    results = calculate_hetkmers(combined_dump, arguments["out_fpath"])
    print(check_run(results))


if __name__ == "__main__":
    main()