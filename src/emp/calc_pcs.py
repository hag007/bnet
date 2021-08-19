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
    parser.add_argument('--phenotype_args', help="phenotype_args", dest='phenotype_args',
                        default=constants.config_json["phenotype_args"])
    args = parser.parse_args()

    calc_pcs(args.dataset_file, args.algo, args.network_file, args.true_solutions_folder, args.additional_args, args.phenotype_args)


def calc_pcs(dataset_file, algo, network_file, true_solutions_folder, additional_args, phenotype_args):

    phenotype_args = json.loads(phenotype_args)
    additional_args = json.loads(additional_args)

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
    phenotype = pd.read_csv(phenotype_args["path to TCGA file"], sep='\t')
    phenotype.index = phenotype.loc[:, constants.config_json["phenotype_index"]]
    phenotype=phenotype.loc[:,phenotype_args["name of variable"]]
    phenotype_values=phenotype_args["value1"]+phenotype_args["value2"]
    phenotype=phenotype[phenotype.isin(phenotype_values)]

    cancer_type=dataset_name.split("_")[1].upper()
    gene_expression_file='/'.join(phenotype_args["path to TCGA file"].split('/')[:-1])+f'/TCGA-{cancer_type}.htseq_fpkm.tsv'
    ge = pd.read_csv(gene_expression_file, sep='\t', index_col=0).transpose()
    gene_names=pd.read_csv(constants.config_json["gene_id_file"],index_col=0).index.values
    ge.columns=[a.split(".")[0] for a in ge.columns]
    ge=ge.reindex(gene_names,axis=1).dropna(axis=1)
    ge=ge.reindex(phenotype.index).dropna(axis=0)


    df_features = pd.DataFrame() # index=ge.index, columns=list(np.arange(n_components)))
    path_to_modules=os.path.join(output_folder,'modules')
    file_names = [f for f in os.listdir(path_to_modules) if os.path.isfile(os.path.join(path_to_modules, f))]
    # iterates over gene files
    print(f'path_to_modules: {path_to_modules}')
    for i, cur_file in enumerate(file_names):

        # print(f'{output_folder}\n{cur_file}\n{constants.config_json["algo"]}_module_genes')

        # print(f'{algo}_module_genes in {cur_file}')

        if f'{algo}_module_genes' in cur_file:

            # print(f'got {algo}_module_genes')

            cur_file_index = os.path.splitext(cur_file)[0].split('_')[-1]

            # reads gene identifiers
            module_genes=pd.read_csv(os.path.join(path_to_modules,cur_file),index_col=0).index.values

            # concatenates PCs of different modules: a SINGLE PC per module
            df_features=pd.concat([df_features,pd.DataFrame(index=ge.index,data=calc_pca(ge.reindex(module_genes, axis=1).dropna(axis=1),n_components=n_components),columns=[cur_file_index])], axis=1)
    print(f'conf: {additional_args}\n# of features: {df_features.shape[1]}')
    df_features.to_csv(output_file, sep='\t')


if __name__ == "__main__":
    main()