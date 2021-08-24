import sys
sys.path.insert(0, '../..')
import os
import pandas as pd


def main(dataset_file, algo, network_file, true_solutions_folder, ts, min_temp, temp_factor, slice_threshold,
         module_threshold, sim_factor, activity_baseline):

    additional_args = [ts, min_temp, temp_factor, slice_threshold, module_threshold, sim_factor,
                             activity_baseline,]

    params_name = "_".join([str(a) for a in additional_args])
    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    network_name = os.path.splitext(os.path.split(network_file)[1])[0]
    output_folder=os.path.join(true_solutions_folder, "{}_{}_{}_{}/modules".format(dataset_name,network_name,algo,params_name))
    output_file=os.path.join(output_folder, "modules_summary.tsv")

    modules = pd.read_csv(output_file, sep='\t')
    active_genes = pd.read_csv(dataset_file, sep='\t')[:len(modules)]['id'].transpose()
    return [[a] for a in active_genes]

