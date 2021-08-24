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
import time

from multiprocessing.pool import Pool
from multiprocessing import Process

from src.emp.generate_solution import generate_solution
from src.emp.calc_pcs import calc_pcs
from src.emp.calc_prediction import calc_prediction

from src.utils.daemon_multiprocessing import MyPool

current_path=os.path.dirname(os.path.realpath(__file__))


def execute_stage(py_function, params):

    # print("about to start script {} with params:\n{}".format(py_function, params))
    py_function(**params)
    # prc=subprocess.Popen("{}../bnet-env/bin/python {} {}".format(constants.dir_path, py_function, params), shell=True,
    #                      stdout=subprocess.PIPE, cwd=current_path)
    # # out = prc.stdout.read()
    # # print(out)
    # while True:
    #     output = prc.stdout.readline()
    #     if output == b'': #  and prc.poll() is not None:
    #         break
    #     if output:
    #         out_str=output.decode("utf-8")
    #         print(out_str, end ='')
    #
    # rc = prc.poll()
    #
    #
    # return rc # 0 # prc.close()



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

    phenotypes_df = pd.read_csv(phenotypes_file, sep='\t')

    algo_param="{}".format(algo)
    true_solutions_folder_param = "{}".format(true_solutions_folder)
    go_folder_param = "{}".format(go_folder)
    pf_param = "{}".format(pf)
    tuning_comb = json.loads(str(tuning_comb_args))
    combs=list(itertools.product(*tuning_comb.values()))
    p=Pool(int(pf))
    params=[]
    for comb in combs:
        for i, row in phenotypes_df.iterrows():
            dataset_file_params = "{}".format(row.loc["path to genes"])
            tcga_path = row.loc["path to TCGA file"]
            var_name = row.loc["name of variable"]
            val1 = row.loc['value1'].replace("'",'\\\"')
            val2 = row.loc['value2'].replace("'",'\\\"')
            # phenotypes_params = '--phenotype_args "{\"path to TCGA file\": \"{}\", \"name of variable\": \"{}\", \"value1\" : [{}], \"value2\" : [{}]}"'.format(\
            #                                     tcga_path, var_name,val1,val2)

            phenotypes_params = json.loads('"{\\\"path to TCGA file\\\": \\\"'+tcga_path+'\\\", \\\"name of variable\\\": \\\"'+var_name+'\\\", \\\"value1\\\" : '+val1+', \\\"value2\\\" : '+val2+'}"')

            prs=[]
            for network_file in network_files:
                network_file_param = network_file

                # pr=Process(target=execute_one_series, args=([additional_args, phenotypes_params, algo_param, comb, dataset_file_params, go_folder_param,
                # network_file_param, processes, true_solutions_folder_param, tuning_comb],))
                # pr.start()
                # prs.append(pr)

                # if len(prs)==1:
                #     for pr in prs:
                #         pr.join()
                #     prs=[]
                #
                params.append([additional_args, phenotypes_params, algo_param, comb, dataset_file_params, go_folder_param,
                network_file_param, processes, true_solutions_folder_param, tuning_comb])


    p.map(execute_one_series,params)


def execute_one_series(args):
    additional_args, phenotypes_params, algo_param, comb, dataset_file_params, go_folder_param, network_file_param, processes, true_solutions_folder_param, tuning_comb = args
    tuning_args = {k: v for k, v in zip(tuning_comb.keys(), comb)}
    additional_args = json.loads(additional_args)
    additional_args.update(tuning_args)
    additional_args['cancer_type']=os.path.basename(dataset_file_params).split('_')[1].upper()
    additional_args['slices_file'] = f'{"/".join(network_file_param.split("/")[:-1])}/{os.path.splitext(os.path.basename(network_file_param))[0]}_louvain_slices.txt'
    additional_args = json.dumps(additional_args)
    additional_args_param = str(additional_args)
    # params_by_process = {
    #     "generate_solution": [dataset_file_params, algo_param, network_file_param, go_folder_param,
    #                           true_solutions_folder_param, additional_args_param],
    #     "calc_pcs": [dataset_file_params, algo_param, network_file_param, go_folder_param,
    #                  true_solutions_folder_param, additional_args_param, phenotypes_params],
    #     "calc_prediction": [dataset_file_params, algo_param, network_file_param, go_folder_param,
    #                         true_solutions_folder_param, additional_args_param, phenotypes_params]}

    function_by_process ={
        "generate_solution": generate_solution,
        "calc_pcs": calc_pcs,
        "calc_prediction": calc_prediction
    }

    params_by_process = {
        "generate_solution": {"dataset_file" : dataset_file_params, "algo" : algo_param, "network_file" : network_file_param , "go_folder" : go_folder_param,
                              "true_solutions_folder": true_solutions_folder_param, "additional_args" : additional_args_param},
        "calc_pcs": {"dataset_file" : dataset_file_params, "algo" : algo_param, "network_file" : network_file_param ,
                              "true_solutions_folder": true_solutions_folder_param, "additional_args" : additional_args_param, "phenotype_args":phenotypes_params},
        "calc_prediction": {"dataset_file" : dataset_file_params, "algo" : algo_param, "network_file" : network_file_param ,
                              "true_solutions_folder": true_solutions_folder_param, "additional_args" : additional_args_param, "phenotype_args":phenotypes_params}}

    dataset_name = os.path.splitext(os.path.split(dataset_file_params)[1])[0]
    network_name = os.path.splitext(os.path.split(network_file_param)[1])[0]
    params_name = "_".join([str(json.loads(additional_args)[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold",
                             "sim_factor", "activity_baseline"]])
    unique_folder_name = "{}_{}_{}_{}".format(dataset_name, network_name, algo_param,
                                              params_name)

    # print(os.path.join(true_solutions_folder_param, unique_folder_name, "report.tsv"))
    if os.path.exists(os.path.join(true_solutions_folder_param, unique_folder_name, "report.tsv")) and constants.config_json['use_cache']=='true':
        print("already exists")
        return

    for cur_process in processes:

        if cur_process not in params_by_process:
            print("unknown process detected: {}. abort...".format(cur_process))
            raise Exception
    try:
        for cur_process in processes:
            execute_stage(function_by_process[cur_process], params_by_process[cur_process])

    except Exception as e:
        print("error in {}: {}".format(cur_process, e))
        raise


if __name__=="__main__":
    main()
