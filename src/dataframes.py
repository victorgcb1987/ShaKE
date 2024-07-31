import pandas as pd

from functools import reduce
from math import log as ln
from math import log2 as log



def merge_kmer_dataframes(dataframes):
    return reduce(lambda df1,df2: pd.merge(df1,df2,on='kmer', how="outer"), dataframes).fillna(0)


def calculate_kmer_count_diversity(dataframe):
    diversity = {"_KmerCount": []}
    col_names = [colname for colname in dataframe.columns if "_kmer_count" in colname and "Specifity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        values = [(float(raw_value)/N) * ln(float(raw_value)/N) if raw_value > 0 else 0 for raw_value in raw_values]
        diversity_value =  -sum(value for value in values if value != 0)
        diversity["_KmerCount"].append(diversity_value)
    dataframe.insert(dataframe.columns.get_loc(col_names[-1])+1, 
                     "Diversity_KmerCount", diversity["_KmerCount"])
    return dataframe


def calculate_kmer_count_specifity(dataframe):
    specifity = {"_KmerCount": []}
    col_names = [colname for colname in dataframe.columns if "_kmer_count" in colname and "Diversity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
        pi = float((1/len(raw_values))) * sum(pijs)
        values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
        si = (1/len(raw_values)) * sum(values)
        specifity["_KmerCount"].append(si)
    dataframe.insert(dataframe.columns.get_loc("Diversity_KmerCount")+1, 
                     "Specifity_KmerCount", specifity["_KmerCount"])
    return dataframe


def calculate_sample_kmer_diversity(dataframe):
    diversity = {"_Sample": []}
    col_names = [colname for colname in dataframe.columns if "kmer" not in colname and "Specifity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        values = [(float(raw_value)/N) * ln(float(raw_value)/N) if raw_value > 0 else 0 for raw_value in raw_values]
        diversity_value =  -sum(value for value in values if value != 0)
        diversity["_Sample"].append(diversity_value)
    dataframe.insert(dataframe.columns.get_loc(col_names[-1])+1, 
                     "Sample_Diversity", diversity["_Sample"])

    return dataframe

    
    
def calculate_sample_kmer_specifity(dataframe):
    specifity = {"_Sample": []}
    col_names = [colname for colname in dataframe.columns if "kmer" not in colname and "Diversity" not in colname]
    for index in dataframe.index:
        raw_values = [dataframe[col_name][index] for col_name in col_names]
        N = sum(raw_values)
        pijs = [abs(float(raw_value/N)) if N > 0 else 0 for raw_value in raw_values]
        pi = float((1/len(raw_values))) * sum(pijs)
        values = [(pij/pi) * log(pij/pi) if pi > 0 and pij > 0 else 0 for pij in pijs]
        si = (1/len(raw_values)) * sum(values)
        specifity["_Sample"].append(si)
    dataframe.insert(dataframe.columns.get_loc("Sample_Diversity")+1, 
                     "Sample_Specifity", specifity["_Sample"])
    return dataframe

def get_kmer_counts_dataframe(filepath, hetkmers=""):
    name = str(filepath.stem).split("_")[0]
    count_colname = "{}_kmer_count".format(name)
    df = pd.read_csv(filepath, delimiter="\t", names=["kmer", "COUNT"])
    if hetkmers:
         for hetkmer in hetkmers:
              df = df.replace(hetkmer[1], hetkmer[0])
         df = df.groupby(["kmer"]).COUNT.sum().reset_index()
    df.columns = ["kmer", count_colname]
    df = df.sort_values(by=["{}_kmer_count".format(name)], ascending=False)
    return df


def rename_dataframe(dataframe):
    kmers = list(dataframe["kmer"])
    ks = ["k{}".format(num[0]+1) for num in enumerate(kmers)]
    index_ = {pair[0]: pair[1] for pair in zip(ks, kmers)}
    dataframe["kmer"] = ks
    return dataframe, index_



