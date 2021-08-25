import sys

sys.path.insert(0, '../../')
# sys.path.insert(0, '../')
# sys.path.insert(0, '../../../')

import os
import subprocess
import src.constants as constants

import pandas as pd
import argparse
import json
import itertools
import numpy as np
import matplotlib.pyplot as plt

from multiprocessing.pool import Pool

current_path = os.path.dirname(os.path.realpath(__file__))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_files', dest='dataset_files', default=constants.config_json["dataset_files"])
    parser.add_argument('--phenotypes_file', dest='phenotypes_file', default=constants.config_json["phenotypes_file"])
    parser.add_argument('--algos', dest='algos', default=constants.config_json["algos"])
    parser.add_argument('--network_file', dest='network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder',
                        default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--metrics_folder', dest='metrics_folder', default=constants.config_json["metrics_folder"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--network_files', dest='network_files', help="network_files", default=constants.config_json["network_files"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args',
                        default=constants.config_json["additional_args"])
    parser.add_argument('--processes', dest="processes", help='generate_solution,calc_pcs,calc_prediction',
                        default=constants.config_json["processes"])
    parser.add_argument('--tuning_comb', dest="tuning_comb", help='',
                        default=constants.config_json["tuning_comb"])
    args = parser.parse_args()
    phenotypes_file = args.phenotypes_file
    algos = args.algos.split(",")
    true_solutions_folder = args.true_solutions_folder
    metrics_folder = args.metrics_folder
    network_files = args.network_files

    network_names = []
    for network_file in network_files:
        network_name = os.path.splitext(os.path.split(network_file)[1])[0]
        network_names.append(network_name)

    metrics=["RF accuracy", "RF AUPR", "RF AUROC","SVM accuracy", "SVM AUPR", "SVM AUROC"]
    null_metrics=["Null SVM AUPR","Null SVM AUROC", "Null RF AUPR","Null RF AUROC"]
    diff_metrics=["SVM AUPR diff","SVM AUROC diff", "RF AUPR diff","RF AUROC diff"]

    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_single_arg_max.tsv")  # suffix_str
    df_algo_max_comb = pd.read_csv(fname, sep='\t', index_col=0)
    df_algo_max_comb = df_algo_max_comb.transpose()
    df_algo_test_comb = df_algo_max_comb.copy()
    df_algo_test_comb.loc[:, :] = None

    for algo in algos:

        suffix_str = f"report_{algo}"
        suffix_str += "_" + network_name
        phenotypes_df = pd.read_csv(phenotypes_file, sep='\t')
        dataset_files=phenotypes_df.loc[:,"path to genes"]




        for metric in metrics:
            if pd.isnull(df_algo_max_comb.loc[algo, metric]):
                continue

            comb_str="_".join(df_algo_max_comb.loc[algo, metric].split("_")[:-1])
            suffix_str += "_" + comb_str

            params = []
            params_null = []
            for dataset_file in dataset_files:
                dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
                suffix_str += "_" + dataset_name

                params.append([dataset_name, network_name, algo, true_solutions_folder]+[comb_str, metric])

                if "AUPR" in metric:
                    params_null.append([dataset_name, network_name, algo, true_solutions_folder] + [comb_str, f"Null {metric}"])
                if "AUROC" in metric:
                    params_null.append([dataset_name, network_name, algo, true_solutions_folder] + [comb_str, f"Null {metric}"])


            results=[fetch_metrics_by_name(prm) for prm in params]
            results=[a for a in results if not a is None]
            df_algo_test_comb.loc[algo, metric]=np.nanmean(results)

            if len(params_null) != 0:
                results = [fetch_metrics_by_name(prm) for prm in params_null]
                results = [a for a in results if not a is None]
                df_algo_test_comb.loc[algo, f"Null {metric}"] = np.mean(results)

            if "AUPR" in metric:
                df_algo_test_comb.loc[algo, f"{metric} diff"]=df_algo_test_comb.loc[algo, f"{metric}"]-df_algo_test_comb.loc[algo, f"Null {metric}"]
            if "AUROC" in metric:
                df_algo_test_comb.loc[algo, f"{metric} diff"]=df_algo_test_comb.loc[algo, f"{metric}"]-df_algo_test_comb.loc[algo, f"Null {metric}"]

    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_test_mean.tsv")  # suffix_str
    df_algo_test_comb.to_csv(fname, sep='\t')



def fetch_metrics_by_name(args):
    dataset_name, network_name, algo, true_solutions_folder,combs_str, metric = args
    # true_solutions_folder, dataset_name, network_name, row ,col, args = args


    result_folder = os.path.join(true_solutions_folder, "{}_{}_{}_{}".format(dataset_name, network_name, algo, combs_str))
    report_file = os.path.join(result_folder, "report.tsv")
    if os.path.exists(report_file):
        print(report_file)
        return pd.read_csv(report_file, sep='\t', index_col=0).fillna(0).iloc[0].loc[metric]
    else:
        return None


def fetch_metrics(args):

    dataset_name, network_name, algo, true_solutions_folder, tuning_args = args

    params_name = "_".join([str(tuning_args[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor",
                             "activity_baseline"]])

    result_folder = os.path.join(true_solutions_folder, "{}_{}_{}_{}".format(dataset_name, network_name, algo, params_name))
    report_file = os.path.join(result_folder, "report.tsv")
    if os.path.exists(report_file):
        df = pd.read_csv(report_file, sep='\t', index_col=0)
    else:
        df=None

    return df


if __name__ == "__main__":
    main()
