import gzip
import os
import struct

from csv import DictReader
from src.utils import get_universe_size


'''BINARY_LENGTH=30(a thousand of millions)!!!!
    Steps (expression):
        1.-Get Normalized data (TPM)
        2.-Round data and convert each value to binary
        3.-Compress data
        4.-Compressed size/original size
        
    Steps (kmer_counts)
        1.-Get universe size
        2.-Calculate difference between universe size and number of kmers in your sample
        3.-Get kmer counts and add 0s equal to this difference
        4.-Convert values to binary (value/total)
        5.-Compress data
        6.-Compressed size/original size'''


def get_universe_size_difference(filepath, universe_size):
    return universe_size - get_universe_size([filepath])


def convert_to_binary(number):
    return format(int(number), '030b')


def create_kmer_binary_file(in_filepath, out_filepath, num_zeros):
    compressed = "{}.gz".format(out_filepath)
    with gzip.open(compressed, 'wb') as compressed_fhand:
        with open(out_filepath, "w") as not_compressed_fhand:
            with open(in_filepath) as in_fhand:
                generator = (convert_to_binary(line.rstrip().split()[1]) for line in in_fhand)
                for gen in generator:
                    compressed_fhand.write(gen.encode()+b"\n")
                    compressed_fhand.flush()
                    not_compressed_fhand.write(gen+"\n")
                    not_compressed_fhand.flush()
                for zero in range(num_zeros):
                    compressed_fhand.write(format(0, '030b').encode()+b"\n")
                    compressed_fhand.flush()
                    not_compressed_fhand.write(format(0, '030b')+"\n")
                    not_compressed_fhand.flush()
    return compressed


def calculate_kolmogorov(filepath_a, filepath_b):
    return float(os.stat(filepath_a).st_size/ os.stat(filepath_b).st_size)


def create_expression_binary_file(in_filepath, units, exclude, out_fpath):
    compressed = "{}.gz".format(out_fpath)
    with gzip.open(compressed, 'wb') as compressed_fhand:
        with open(out_fpath, "w") as not_compressed_fhand:
            with open(in_filepath) as fhand:
                generator = (convert_to_binary(round(float(row[units]), 3)*1000) for row in DictReader(fhand, delimiter="\t") if row["Reference"] not in exclude)
                for gen in generator:
                    compressed_fhand.write(gen.encode()+b"\n")
                    compressed_fhand.flush()
                    not_compressed_fhand.write(gen+"\n")
                    not_compressed_fhand.flush()
    return compressed

def calculate_kolmogorov_estimator(filepath, universe_size, estimators, group=None, 
                                   sub=None, name=None, kind=None, units="TPM"):
    binary = "{}.binary".format(str(filepath))
    if kind != "expression":
        num_zeros = get_universe_size_difference(filepath, universe_size)
        compressed_file = create_kmer_binary_file(filepath, binary, num_zeros)
    else:
        compressed_file = create_expression_binary_file(filepath, units, [], binary)
    kolmo = calculate_kolmogorov(compressed_file, binary)
    estimators[group][sub][name]["kolmogorov"] = kolmo


    