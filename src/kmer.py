import pandas as pd

from io import StringIO
from math import log as ln
from src.utils import reformat_lines


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

def group_kmers_counts(filepaths, hetkmers={}):
    kmers_counts = {"header": []}
    kmer_universe = []
    sample_values = []
    for filepath in filepaths:
          kmers_counts["header"].append(filepath.stem.split(".")[0])
          with open(filepath) as filehand:
               values = {val[0]: val[1] for val in {(hetkmers.get(count[0], count[0]), int(count[1])) for count in (count.split() for count in filehand.read().split("\n")) if len(count) == 2}}
               sample_values.append(values)
               kmer_universe += list(values.keys())
    for kmer in set(kmer_universe):
          kmers_counts[kmer] = [val.get(kmer, 0) for val in sample_values]
    return kmers_counts
          

def index_kmers(kmer_counts):
     n = 0
     index = []
     kmers = list(kmer_counts.keys())
     for kmer in kmers:
          if kmer == "header":
               continue
          n += 1
          indx = "K{}".format(n)
          index.append((indx, kmer))
          kmer_counts[indx] = kmer_counts.pop(kmer)
     return index, kmer_counts 


# def calculate_sample_estimators(kmer_counts):
#      sample_diversity = {sample: 0 for sample in kmer_counts["header"]}
#      for pos, samplename in enumerate(kmer_counts["header"]):
#           raw_values = [kmer_counts[kmer][pos] for kmer in kmer_counts.keys() if kmer != "header"]
#           N = sum(raw_values)
#           values = [(float(value)/N) * ln(float(value)/N) if value > 0 else 0 for value in raw_values]
#           diversity_value =  -sum(value for value in values if value != 0)
#           sample_diversity[samplename] = diversity_value
#      sample_specifity = {sample: 0 for sample in kmer_counts["header"]}
#      for pos, samplename in enumerate(kmer_counts["header"]):
#           raw_values = [kmer_counts[kmer][pos] for kmer in kmer_counts.keys() if kmer != "header"]
#           N = sum(raw_values)
#           pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
#           pi = float((1/len(raw_values))) * sum(pijs)
#           values = [(pij/pi) * ln(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
#           si = (1/len(raw_values)) * sum(values)
#           sample_specifity[samplename] = si
#      return sample_diversity, sample_specifity
     

def calculate_kmer_estimators(kmer_counts):
     kmer_diversity = {kmer: 0 for kmer in kmer_counts if kmer != "header"}
     for kmer, counts in kmer_counts.items():
          if kmer == "header":
               continue
          raw_values = [count for count in counts]
          N = sum(raw_values)
          values = [(float(value)/N) * ln(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_value =  -sum(value for value in values if value != 0)
          kmer_diversity[kmer] = diversity_value
     kmer_specifity = {kmer: 0 for kmer in kmer_counts if kmer != "header"}
     for kmer, counts in kmer_counts.items():
          if kmer == "header":
               continue
          raw_values = [count for count in counts]
          N = sum(raw_values)
          pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
          pi = float((1/len(raw_values))) * sum(pijs)
          values = [(pij/pi) * ln(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
          si = (1/len(raw_values)) * sum(values)
          kmer_specifity[kmer] = si
     return kmer_diversity, kmer_specifity

def calculate_sample_estimators(filepath, universe_size, estimators):
     with open(filepath) as fhand:
          raw_values = [int(line.rstrip().split()[1]) for line in fhand if line]
          N = sum(raw_values)
          #diversity
          values = [(float(value)/N) * ln(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_value =  -sum(value for value in values if value != 0)
          #especifity
          pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
          pi = float((1/len(raw_values))) * sum(pijs)
          values = [(pij/pi) * ln(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
          especifity = (1/universe_size) * sum(values)
          estimators[filepath.stem] = {"diversity": diversity_value, "especifity": especifity}
