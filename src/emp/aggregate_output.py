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

from multiprocessing.pool import Pool

current_path = os.path.dirname(os.path.realpath(__file__))


def main():
    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_files', dest='dataset_files', default=constants.config_json["dataset_files"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder',
                        default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--metrics_folder', dest='metrics_folder', default=constants.config_json["metrics_folder"])
    parser.add_argument('--pf', dest='pf', help="parallelization factor", default=constants.config_json["pf"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args',
                        default=constants.config_json["additional_args"])
    parser.add_argument('--processes', dest="processes", help='generate_solution,calc_pcs,calc_prediction',
                        default=constants.config_json["processes"])
    parser.add_argument('--tuning_comb', dest="tuning_comb", help='',
                        default=constants.config_json["tuning_comb"])
    args = parser.parse_args()
    dataset_files = args.dataset_files
    algo = args.algo
    true_solutions_folder = args.true_solutions_folder
    tuning_comb_args = args.tuning_comb
    metrics_folder = args.metrics_folder

    tuning_comb = json.loads(str(tuning_comb_args))
    combs = list(itertools.product(*tuning_comb.values()))
    p = Pool(30)
    params = []
    combs_str = f"report_{algo}"
    for dataset_file in dataset_files:
        dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
        combs_str += "_" + dataset_name
        for comb in combs:
            params.append([dataset_name, algo, true_solutions_folder, comb, tuning_comb])

    combs_str += "_" + "_".join([str(tuning_comb[a]) for a in \
                                                                ["ts", "min_temp", "temp_factor", "slice_threshold",
                                                                 "module_threshold", "sim_factor",
                                                                 "activity_baseline"]])
    results = p.map(fetch_metrics, params)

    df = pd.DataFrame()
    for result in results:
        df = df.append(result)

    if not os.path.isdir(metrics_folder):
        os.makedirs(metrics_folder)
    fname = os.path.join(metrics_folder, f"{combs_str}.csv")
    df.to_csv(fname)


def fetch_metrics(args):
    dataset_name, algo, true_solutions_folder, comb, tuning_comb = args
    tuning_args = {k: v for k, v in zip(tuning_comb.keys(), comb)}
    params_name = "_".join([str(tuning_args[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor",
                             "activity_baseline"]])

    result_folder = os.path.join(true_solutions_folder, "{}_{}_{}".format(dataset_name, algo, params_name))
    report_file = os.path.join(result_folder, "report.tsv")
    df = pd.read_csv(report_file,index_col=0)
    return df


if __name__ == "__main__":
    main()
