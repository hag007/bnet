import sys
sys.path.insert(0, '../..')
import os
import pandas as pd


def main(dataset_file, compare_folder):
    output_folder = os.path.join(compare_folder, "modules")
    output_file=os.path.join(output_folder, "modules_summary.tsv")

    modules = pd.read_csv(output_file, sep='\t')
    active_genes = pd.read_csv(dataset_file, sep='\t')[:len(modules)]['id'].transpose()
    return [[a] for a in active_genes]

