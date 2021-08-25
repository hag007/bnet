import sys

sys.path.insert(0, '../../')

import os
import src.constants as constants

import pandas as pd
import argparse
import json
import itertools

from run_sequentially import run_sequentially

current_path = os.path.dirname(os.path.realpath(__file__))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_files', dest='dataset_files', default=constants.config_json["dataset_files"])
    parser.add_argument('--phenotypes_file', dest='phenotypes_file', default=constants.config_json["phenotypes_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_files', dest='network_files', default=constants.config_json["network_files"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default=constants.config_json["additional_args"])
    parser.add_argument('--processes', dest="processes", help='generate_solution,calc_pcs,calc_prediction', default=constants.config_json["processes"])
    parser.add_argument('--tuning_comb', dest="tuning_comb", help='', default=constants.config_json["tuning_comb"])
    parser.add_argument('--algos', dest='algos', default=constants.config_json["algos"])
    parser.add_argument('--metrics_folder', dest='metrics_folder', default=constants.config_json["metrics_folder"])
    args = parser.parse_args()

    dataset_files=args.dataset_files
    phenotypes_file=args.phenotypes_file
    algo=args.algo
    network_files=args.network_files
    network_name = os.path.splitext(os.path.split(network_files[0])[-1])[0]
    go_folder=args.go_folder
    true_solutions_folder=args.true_solutions_folder
    pf=args.pf
    additional_args = args.additional_args
    tuning_comb_args = args.tuning_comb
    processes = args.processes.split(",")
    algos = args.algos.split(",")
    metrics_folder = args.metrics_folder



    fields=["N features", "RF accuracy", "RF AUPR", "RF AUROC", "Null RF AUPR", "Null RF AUROC", "SVM accuracy","SVM AUPR","SVM AUROC","Null SVM AUPR","Null SVM AUROC","SVM AUPR diff","SVM AUROC diff", "number of datasets"]

    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_single_arg_max.tsv")  # combs_str

    df_single_arg_max=pd.to_csv(fname, sep='\t', index_col=0).iloc[:,0]

    for algo in algos:

        for comapre_algo,comb in df_single_arg_max.iterrows():
            comb_json={}
            comb_arr=comb.split("_")
            comb_json["ts"]=comb_arr[0]
            comb_json["min_temp"] = comb_arr[1]
            comb_json["temp_factor"] = comb_arr[2]
            comb_json["slice_threshold"] = comb_arr[3]
            comb_json["module_threshold"] = comb_arr[4]
            comb_json["sim_factor"] = comb_arr[5]
            comb_json["activity_baseline"] = comb_arr[6]

            run_sequentially(algo,network_files,constants.config_json["phenotypes_file"],constants.config_json["processes"],constants.config_json["go_folder"],
                             constants.config_json["true_solutions_folder"],constants.config_json["pf"],constants.config_json["additiona_args"],json.dumps(comb_json))



if __name__ == "__main__":
    main()
