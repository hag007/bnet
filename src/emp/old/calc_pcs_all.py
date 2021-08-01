import sys
sys.path.insert(0, '../../..')
import json
import argparse
import os
import src.constants as constants
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
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
    additional_args = args.additional_args

    dataset_name=os.path.splitext(os.path.split(dataset_file)[1])[0]
    output_folder=os.path.join(true_solutions_folder, "{}_{}".format(dataset_name,algo))
    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    n_components=4

    phenotype = pd.read_csv(constants.config_json["phenotype_file"], sep='\t')
    phenotype.index = phenotype.loc[:, constants.config_json["phenotype_index"]]
    phenotype=phenotype.loc[:,constants.config_json["phenotype_field"]]
    phenotype=phenotype[phenotype.isin(constants.config_json["phenotype_values"])]

    ge = pd.read_csv(constants.config_json["gene_expression_file"], sep='\t', index_col=0).transpose()
    gene_names=pd.read_csv(constants.config_json["gene_id_file"],index_col=0).index.values
    ge.columns=[a.split(".")[0] for a in ge.columns]
    ge=ge.loc[:,gene_names].dropna(axis=1)
    ge=ge.reindex(phenotype.index).dropna(axis=0)


    path_to_modules=f'{constants.config_json["true_solutions_folder"]}/{os.path.splitext(os.path.basename(constants.config_json["dataset_file"]))[0]}_{constants.config_json["algo"]}/report/'
    df_features=pd.DataFrame(index=ge.index,data=calc_pca(ge,n_components=n_components),columns=list(np.arange(4)))
    df_features.to_csv(constants.config_json["pcs_file"], sep='\t')


if __name__ == "__main__":
    main()