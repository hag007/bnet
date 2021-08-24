import sys
sys.path.insert(0, '../..')
import json
import argparse
import os
import src.constants as constants
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from src.emp.calc_prediction import main as calc_prediction_main

# calc first PC
def calc_pca(X, n_components=1):
    pca = PCA(n_components=n_components)
    X=pca.fit_transform(X.values)
    # X=np.mean(X.values, axis=1).transpose()
    return X

def main(dataset_file, algo, network_file, true_solutions_folder, ts, min_temp, temp_factor, slice_threshold,
         module_threshold, sim_factor, activity_baseline):

    additional_args = [ts, min_temp, temp_factor, slice_threshold, module_threshold, sim_factor,
                       activity_baseline, ]

    params_name = "_".join([str(a) for a in additional_args])
    dataset_name = os.path.splitext(os.path.split(dataset_file)[1])[0]
    network_name = os.path.splitext(os.path.split(network_file)[1])[0]
    orig_true_solutions = "/home/gaga/hagailevi/omics/output/true_solutions"
    orig_output_folder = os.path.join(orig_true_solutions,
                                 "{}_{}_{}_{}/modules".format(dataset_name, network_name, algo, params_name))
    orig_output_file = os.path.join(orig_output_folder, "modules_summary.tsv")

    output_folder = os.path.join(true_solutions_folder,
                                 "{}_{}_{}_{}/modules".format(dataset_name, network_name, algo, params_name))

    modules = pd.read_csv(orig_output_file, sep='\t')
    num_modules = len(modules)
    num_of_genes = sum(modules['#_genes'])
    active_genes = pd.read_csv(dataset_file, sep='\t')[:num_of_genes].transpose()
    active_genes = active_genes.transpose()['id']

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    cancer_type=dataset_name.split("_")[1].upper()
    path_to_TCGA_file = f'/specific/netapp5/gaga/hagailevi/omics/GDC-TCGA/{cancer_type}/tcga_data/TCGA-{cancer_type}.GDC_phenotype.tsv'
    gene_expression_file='/'.join(path_to_TCGA_file.split('/')[:-1])+f'/TCGA-{cancer_type}.htseq_fpkm.tsv'
    ge = pd.read_csv(gene_expression_file, sep='\t', index_col=0).transpose()
    gene_names=pd.read_csv(constants.config_json["gene_id_file"],index_col=0).index.values
    ge.columns=[a.split(".")[0] for a in ge.columns]
    ge=ge.reindex(gene_names,axis=1).dropna(axis=1)
    ge_relevant = ge[active_genes]
    pca = calc_pca(ge_relevant, num_modules)
    output_file = os.path.join(output_folder, "pcs.tsv")
    pd.DataFrame(pca).to_csv(output_file, sep='\t', index=False)
    os.chmod(output_file, 0o777)

