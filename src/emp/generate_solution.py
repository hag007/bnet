import sys
sys.path.insert(0, '../..')
# import pandas as pd
import json
import argparse
import os
from src.runners.run_algo import run_algo
import src.constants as constants
from src.utils.go import init_state


def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--compare', dest='compare', default=constants.config_json["compare_algo"])
    parser.add_argument('--network_file', dest='network_file', help='/path/to/network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default=constants.config_json["additional_args"])
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    compare_algo = args.compare
    network_file = args.network_file
    go_folder = args.go_folder
    true_solutions_folder = args.true_solutions_folder
    additional_args = args.additional_args

    generate_solution(dataset_file, algo, go_folder, network_file, true_solutions_folder, compare_algo, additional_args)


def generate_solution(dataset_file, algo, go_folder, network_file, true_solutions_folder, compare_algo, additional_args):

    print(f'start ({(dataset_file, algo, go_folder, network_file, true_solutions_folder ,additional_args)})')
    # init_state(go_folder)
    additional_args_json = json.loads(additional_args)
    params_name = "_".join([str(additional_args_json[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor",
                             "activity_baseline"]])
    # read files
    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
    network_name = os.path.splitext(os.path.split(network_file)[1])[0]
    output_folder = os.path.join(true_solutions_folder,
                                 "{}_{}_{}_{}".format(dataset_name, network_name, algo, params_name))
    if compare_algo is not None:
        compare_folder = os.path.join(true_solutions_folder, "{}_{}_{}_{}".format(dataset_name,network_name,compare_algo,params_name))
        additional_args_json['compare'] = compare_folder

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    run_algo(dataset_file, algo, network_file, go_folder, os.path.join(output_folder) , **additional_args_json)
