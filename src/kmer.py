from math import log2 as log

from src.utils import get_kmer_value


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

def calculate_kmer_estimators(kmer_counts):
     kmer_diversity = {kmer: 0 for kmer in kmer_counts if kmer != "header"}
     for kmer, counts in kmer_counts.items():
          if kmer == "header":
               continue
          raw_values = [count for count in counts]
          N = sum(raw_values)
          values = [(float(value)/N) * log(float(value)/N) if value > 0 else 0 for value in raw_values]
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
          values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
          si = (1/len(raw_values)) * sum(values)
          kmer_specifity[kmer] = si
     return kmer_diversity, kmer_specifity


def calculate_sample_estimators(filepath, universe_size, estimators, group=None, 
                                sub=None, name=None, file=None, pipe=False, kind=None, binary=False):
     with open(filepath) as fhand:
          raw_values = [int(line.rstrip().split()[1]) for line in fhand if line]
          if binary:
               raw_values = [1 if float(value) >= 1 else 0 for value in raw_values]
          N = sum(raw_values)
          #diversity
          values = [(float(value)/N) * log(float(value)/N) if value > 0 else 0 for value in raw_values]
          diversity_value =  -sum(value for value in values if value != 0)
          #especifity
          #pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
          #pi = float((1/len(raw_values))) * sum(pijs)
          #values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
          specifity = log(universe_size) - diversity_value
          results = {"diversity": diversity_value, "specifity": specifity, 
                     "universe_size": universe_size, "sub": sub,
                     "name": name, "kind": kind, "file": file}
          if pipe:
               if group not in estimators:
                    estimators[group] = {}
               if sub not in estimators[group]:
                    estimators[group][sub] = {name: results}
               else:
                    estimators[group][sub][name] = results
          else:
               estimators[filepath.stem] = {"diversity": diversity_value, "specifity": specifity}


def calculate_kmer_estimators(filepaths, universe_size , kmer):
     raw_values = [get_kmer_value(filepath, kmer) for filepath in filepaths]
     N = sum(raw_values)
     values = [(float(value)/N) * log(float(value)/N) if value > 0 else 0 for value in raw_values]
     diversity_value =  -sum(value for value in values if value != 0)
     specifity = log(universe_size) - diversity_value
     return diversity_value, specifity 