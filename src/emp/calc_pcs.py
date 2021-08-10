import sys
sys.path.insert(0, '../..')
import json
import argparse
import os
import src.constants as constants
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

# calc first PC
def calc_pca(X, n_components=1):
    pca = PCA(n_components=n_components)
    X=pca.fit_transform(X.values)
    # X=np.mean(X.values, axis=1).transpose()
    return X

def main():

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--dataset_file', dest='dataset_file', help='/path/to/dataset_file', default=constants.config_json["dataset_file"])
    parser.add_argument('--algo', dest='algo', default=constants.config_json["algo"])
    parser.add_argument('--network_file', dest='network_file', help='/path/to/network_file', default=constants.config_json["network_file"])
    parser.add_argument('--go_folder', dest='go_folder', default=constants.config_json["go_folder"])
    parser.add_argument('--true_solutions_folder', dest='true_solutions_folder', default=constants.config_json["true_solutions_folder"])
    parser.add_argument('--additional_args', help="additional_args", dest='additional_args', default=constants.config_json["additional_args"])
    args = parser.parse_args()

    dataset_file=args.dataset_file
    algo=args.algo
    network_file = args.network_file
    go_folder = args.go_folder
    true_solutions_folder = args.true_solutions_folder

    additional_args = json.loads(args.additional_args)

    params_name = "_".join([str(additional_args[a]) for a in \
                            ["ts", "min_temp", "temp_factor", "slice_threshold", "module_threshold", "sim_factor",
                             "activity_baseline"]])

    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    network_name = os.path.splitext(os.path.split(network_file)[1])[0]
    output_folder=os.path.join(true_solutions_folder, "{}_{}_{}_{}".format(dataset_name,network_name,algo,params_name))
    output_file=os.path.join(output_folder, "pcs.tsv")
    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    n_components=1

    # read datasets
    phenotype = pd.read_csv(constants.config_json["phenotype_file"], sep='\t')
    phenotype.index = phenotype.loc[:, constants.config_json["phenotype_index"]]
    phenotype=phenotype.loc[:,constants.config_json["phenotype_field"]]
    phenotype=phenotype[phenotype.isin(constants.config_json["phenotype_values"])]

    ge = pd.read_csv(constants.config_json["gene_expression_file"], sep='\t', index_col=0).transpose()
    gene_names=pd.read_csv(constants.config_json["gene_id_file"],index_col=0).index.values
    ge.columns=[a.split(".")[0] for a in ge.columns]
    ge=ge.loc[:,gene_names].dropna(axis=1)
    ge=ge.reindex(phenotype.index).dropna(axis=0)


    df_features = pd.DataFrame() # index=ge.index, columns=list(np.arange(n_components)))
    path_to_modules=os.path.join(output_folder,'modules')
    file_names = [f for f in os.listdir(path_to_modules) if os.path.isfile(os.path.join(path_to_modules, f))]
    # iterates over gene files
    for i, cur_file in enumerate(file_names):
        if f'{constants.config_json["algo"]}_module_genes' in cur_file:
            cur_file_index = os.path.splitext(cur_file)[0].split('_')[-1]

            # reads gene identifiers
            module_genes=pd.read_csv(os.path.join(path_to_modules,cur_file),index_col=0).index.values

            # concatenates PCs of different modules: a SINGLE PC per module
            df_features=pd.concat([df_features,pd.DataFrame(index=ge.index,data=calc_pca(ge.loc[:,module_genes].dropna(axis=1),n_components=n_components),columns=[cur_file_index])], axis=1)
    df_features.to_csv(output_file, sep='\t')


if __name__ == "__main__":
    main()