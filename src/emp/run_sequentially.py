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

current_path=os.path.dirname(os.path.realpath(__file__))


def execute_stage(py_script, params):

    params = " ".join(params)
    print("about to start script {} with params:\n{}".format(py_script, params))
    prc=subprocess.Popen("{}../bnet-env/bin/python {} {}".format(constants.dir_path, py_script, params), shell=True,
                           stdout=subprocess.PIPE, cwd=current_path)
    # out = prc.stdout.read()
    # print(out)
    while True:
        output = prc.stdout.readline()
        if output == b'': #  and prc.poll() is not None:
            break
        if output:
            out_str=output.decode("utf-8")
            print(out_str, end ='')

    rc = prc.poll()


    return rc # 0 # prc.close()



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
    parser.add_argument('--tuning_comb', dest="tuning_comb", help='',
                        default=constants.config_json["tuning_comb"])
    args = parser.parse_args()

    dataset_files=args.dataset_files
    phenotypes_file=args.phenotypes_file
    algo=args.algo
    network_files=args.network_files
    go_folder=args.go_folder
    true_solutions_folder=args.true_solutions_folder
    pf=args.pf
    additional_args = args.additional_args
    tuning_comb_args = args.tuning_comb
    processes = args.processes.split(",")

    phenotypes_df = pd.read_csv(phenotypes_file)


    algo_param="--algo {}".format(algo)
    true_solutions_folder_param = "--true_solutions_folder {}".format(true_solutions_folder)
    go_folder_param = "--go_folder {}".format(go_folder)
    pf_param = "--pf {}".format(pf)
    tuning_comb = json.loads(str(tuning_comb_args))
    combs=list(itertools.product(*tuning_comb.values()))
    p=Pool(1)
    params=[]
    for i, row in phenotypes_df.iterrows():
        dataset_file_params = "--dataset_file {}".format(row.loc["path to genes"])
        tcga_path = row.loc["path to TCGA file"]
        var_name = row.loc["name of variable"]
        val1 = row.loc['value1'].replace("'",'\\\"')
        val2 = row.loc['value2'].replace("'",'\\\"')
        # phenotypes_params = '--phenotype_args "{\"path to TCGA file\": \"{}\", \"name of variable\": \"{}\", \"value1\" : [{}], \"value2\" : [{}]}"'.format(\
        #                                     tcga_path, var_name,val1,val2)

        phenotypes_params = '--phenotype_args "{\\\"path to TCGA file\\\": \\\"'+tcga_path+'\\\", \\\"name of variable\\\": \\\"'+var_name+'\\\", \\\"value1\\\" : '+val1+', \\\"value2\\\" : '+val2+'}"'


        for network_file in network_files:
            network_file_param = "--network_file {}".format(network_file)
            for comb in combs:
                params.append([additional_args, phenotypes_params, algo_param, comb, dataset_file_params, go_folder_param,
                network_file_param, processes, true_solutions_folder_param, tuning_comb])
    p.map(execute_one_series,params)


def execute_one_series(args):
    additional_args, phenotypes_params, algo_param, comb, dataset_file_params, go_folder_param, network_file_param, processes, true_solutions_folder_param, tuning_comb = args
    tuning_args = {k: v for k, v in zip(tuning_comb.keys(), comb)}
    additional_args = json.loads(additional_args)
    additional_args.update(tuning_args)
    print(additional_args)
    additional_args = json.dumps(additional_args)
    additional_args_param = "--additional_args {}".format(json.dumps(str(additional_args)))
    params_by_processes = {
        "generate_solution": [dataset_file_params, algo_param, network_file_param, go_folder_param,
                              true_solutions_folder_param, additional_args_param],
        "calc_pcs": [dataset_file_params, algo_param, network_file_param, go_folder_param,
                     true_solutions_folder_param, additional_args_param],
        "calc_prediction": [dataset_file_params, algo_param, network_file_param, go_folder_param,
                            true_solutions_folder_param, additional_args_param, phenotypes_params]}
    for cur_process in processes:
        if cur_process not in params_by_processes:
            print("unknown process detected: {}. abort...".format(cur_process))
            raise Exception
    try:
        for cur_process in processes:
            execute_stage(cur_process + ".py", params_by_processes[cur_process])

    except Exception as e:
        print("error in {}: {}".format(cur_process, e))
        raise


if __name__=="__main__":
    main()
