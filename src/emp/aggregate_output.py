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
    dataset_files = args.dataset_files
    phenotypes_file = args.phenotypes_file
    algos = args.algos.split(",")
    true_solutions_folder = args.true_solutions_folder
    tuning_comb_args = args.tuning_comb
    metrics_folder = args.metrics_folder
    network_files = args.network_files

    tuning_comb = json.loads(str(tuning_comb_args))
    combs = list(itertools.product(*tuning_comb.values()))

    fields=["N features", "RF accuracy", "RF AUPR", "RF AUROC", "Null RF AUPR", "Null RF AUROC", "SVM accuracy","SVM AUPR","SVM AUROC","Null SVM AUPR","Null SVM AUROC","SVM AUPR diff","SVM AUROC diff", "number of datasets"]
    metric_fields = ["RF accuracy", "RF AUPR", "RF AUROC", "SVM accuracy", "SVM AUPR", "SVM AUROC"]
    df_algos_mean=pd.DataFrame(columns=fields)
    df_algos_median = pd.DataFrame(
        columns=fields)
    df_algos_max = pd.DataFrame(
        columns=fields)
    df_algo_max_comb = pd.DataFrame()
    df_ranking= pd.DataFrame()
    for algo in algos:

        # p = Pool(5)
        params = []
        combs_str = f"report_{algo}"

        phenotypes_df = pd.read_csv(phenotypes_file, sep='\t')

        dataset_files=phenotypes_df.loc[:,"path to genes"]

        df_means=pd.DataFrame(columns=fields)
        df_median = pd.DataFrame(
            columns=fields)
        df_max = pd.DataFrame(columns=fields)
        for comb in combs:

            cur_tuning_json = {k: v for k, v in zip(tuning_comb.keys(), comb)}
            params=[]

            network_names=[]
            for network_file in network_files:
                network_name = os.path.splitext(os.path.split(network_file)[-1])[0]
                combs_str += "_" + network_name
                network_names.append(network_name)

                for dataset_file in dataset_files:
                    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
                    combs_str += "_" + dataset_name

                    params.append([dataset_name, network_name, algo, true_solutions_folder, cur_tuning_json])

            cur_tuning_json_str="_".join([str(cur_tuning_json[a]) for a in \
                                                                        ["ts", "min_temp", "temp_factor", "slice_threshold",
                                                                         "module_threshold", "sim_factor",
                                                                         "activity_baseline"]])
            network_names_str = "_".join(network_names)
            combs_str += "_"+cur_tuning_json_str
            # results = p.map(, params)
            results=[fetch_metrics(prm) for prm in params]
            results=[a for a in results  if not a is None]
            df = pd.DataFrame()
            for result in results:
                df = pd.concat([df, result])

            if not os.path.isdir(metrics_folder):
                os.makedirs(metrics_folder)
            fname = os.path.join(metrics_folder, f"agg_report_{algo}_{cur_tuning_json_str}_{network_names_str}.tsv") # combs_str
            print(fname)
            df.to_csv(fname, sep='\t')
            cur_index=f'{cur_tuning_json_str}_{network_names_str}'
            if df_means is None:
                df_means=pd.DataFrame(index=[], data=[df.dropna(axis=0).mean(axis=0)])
            if df_median is None:
                df_median=pd.DataFrame(index=[], data=[df.dropna(axis=0).median(axis=0)])
            else:
                if not df.empty:
                    df.loc[:, "RF AUPR diff"] = df.loc[:, "RF AUPR"] - df.loc[:, "Null RF AUPR"]
                    df.loc[:, "RF AUROC diff"] = df.loc[:, "RF AUROC"] - df.loc[:, "Null RF AUROC"]
                    df.loc[:, "SVM AUPR diff"] = df.loc[:, "SVM AUPR"] - df.loc[:, "Null SVM AUPR"]
                    df.loc[:, "SVM AUROC diff"] = df.loc[:, "SVM AUROC"] - df.loc[:, "Null SVM AUROC"]
                df_means.loc[cur_index]=df.dropna(axis=0).mean(axis=0)
                df_means.loc[cur_index, "number of datasets"] = df.dropna(axis=0).shape[0]
                df_median.loc[cur_index] = df.dropna(axis=0).mean(axis=0)
                df_median.loc[cur_index, "number of datasets"] = df.dropna(axis=0).shape[0]
                df_max.loc[cur_index] = df.dropna(axis=0).max(axis=0)
                df_max.loc[cur_index, "number of datasets"] = df.dropna(axis=0).shape[0]

            # for col in df.columns:
            #     if not col.startswith("Null"):
            #         plt.hist(df[col], bins='auto')  # arguments are passed to np.histogram
            #         plt.title(f"{algo} - {network_names_str} - {col} - {cur_tuning_json_str}\n{np.std(df[col].values)}")
            #
            #         plt.savefig(
            #             f"/home/gaga/hagailevi/omics/output/histograms/{algo} - {network_names_str} - {col} - {cur_tuning_json_str}.png")
            #         plt.close('all')


        fname = os.path.join(metrics_folder, f"agg_report_{algo}_{network_names_str}.tsv")  # combs_str
        df_means=df_means.dropna(axis=0)
        df_means.to_csv(fname, sep='\t')
        df_algos_mean.loc[algo]=df_means.mean(axis=0)
        df_algos_mean.loc[algo,"number of combinations"] = df_means.shape[0]
        df_algos_median.loc[algo] = df_means.mean(axis=0)
        df_algos_median.loc[algo, "number of combinations"] = df_means.shape[0]
        df_algos_max.loc[algo] = df_means.max(axis=0)
        # cur_algo_ranking=df_means.loc[:,["RF accuracy", "RF AUPR", "RF AUROC", "SVM accuracy","SVM AUPR","SVM AUROC"]].rank(axis=0).applymap(lambda a: 1 if a == df_means.shape[0] else 0.5 if a == df_means.shape[0]-1 else 0).sum(axis=1)
        cur_algo_ranking = df_means.loc[:,
                           ["RF accuracy", "RF AUPR", "RF AUROC", "SVM accuracy", "SVM AUPR", "SVM AUROC"]].apply(lambda a:  a.max()-a  ,axis=0).applymap(lambda a: 1-min(a, 0.01)*100).sum(
            axis=1)
        df_ranking=pd.concat([df_ranking, cur_algo_ranking.to_frame().rename(columns={0:algo})], axis=1)

        df_algos_max.loc[algo, "number of combinations"] = df_means.shape[0]

        df_max_comb=pd.isnull(df_means[df_means == df_means.max(axis=0)]).applymap(lambda a: not a)
        df_max_comb=df_max_comb.drop(columns=['N features', 'Null RF AUROC', 'Null RF AUPR', 'Null SVM AUROC', 'Null SVM AUPR', 'number of datasets'])
        df_max_comb=df_max_comb[df_max_comb.apply(lambda a: any(a), axis=1)]
        df_max_comb=df_max_comb.apply(lambda a: list(df_max_comb.index[a]))
        if not df_max_comb.empty:
            df_algo_max_comb.loc[:,algo]=pd.DataFrame(df_max_comb).iloc[0]

    # df_ranking.apply(lambda a:  (a.max()==a).astype(int)  ,axis=1)




    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_single_arg_max.tsv")  # combs_str
    df_single_arg_max=pd.concat([df_ranking.apply(lambda b: b.index.values[b.argmax()], axis=0).to_frame(name=a) for a in metric_fields], axis =1).transpose()
    print(fname)
    df_single_arg_max.to_csv(fname, sep='\t')



    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_mean.tsv")  # combs_str
    df_algos_mean.to_csv(fname, sep='\t')
    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_median.tsv")  # combs_str
    df_algos_median.to_csv(fname, sep='\t')
    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_max.tsv")  # combs_str
    df_algos_max.to_csv(fname, sep='\t')
    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_arg_max.tsv")  # combs_str
    df_algo_max_comb.to_csv(fname, sep='\t')
    fname = os.path.join(metrics_folder, f"agg_report_{network_name}_{'_'.join(algos)}_ranking.tsv")  # combs_str
    df_ranking.to_csv(fname, sep='\t')



def fetch_metrics(args):

    dataset_name, network_name, algo, true_solutions_folder, tuning_args = args

    params_name = "_".join([str(tuning_args[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor",\
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

