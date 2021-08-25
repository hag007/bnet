import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../..')
import os
import pandas as pd

from src import constants
from src.implementations.top_sdg import main as top_sdg_main
from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.network import get_network_genes
from src.utils.go_similarity import init_go_metadata

from src.runners.abstract_runner import AbstractRunner


class TopSDGGenesRunner(AbstractRunner):
    def __init__(self):
        super().__init__(f"top_SDG_genes")

    def extract_modules_and_bg(self, bg_genes, dest_algo_dir):
        results = open(os.path.join(dest_algo_dir, "modules.txt")).readlines()
        modules = [[] for x in range(max([int(x.strip().split(" =")[1]) for x in results[1:]]) + 1)]
        for x in results[1:]:
            if int(x.strip().split(" =")[1]) != -1:
                modules[int(x.strip().split(" =")[1])].append(x.strip().split(" =")[0])
            else:
                modules.append([x.strip().split(" =")[0]])
        modules = filter(lambda x: len(x) > 3, modules)
        all_bg_genes = [bg_genes for x in modules]
        print("extracted {} modules".format(len(modules)))
        return modules, all_bg_genes

    def init_params(self, dataset_file_name, network_file_name, output_folder):


        df_scores=pd.read_csv(dataset_file_name, sep='\t', index_col=0)
        sig_genes=df_scores['qval'][df_scores['qval']<0.05].index
        active_genes_file=os.path.join(output_folder, "active_genes_file.txt")
        open(active_genes_file, "w+").write("\n".join([x for x in sig_genes if len(ensembl2entrez_convertor([x]))>0 ]))
        bg_genes=get_network_genes(network_file_name)
        return active_genes_file, bg_genes


    def run(self, dataset_file_name, network_file_name, output_folder, **kwargs):
        print("run top_sdg_genes runner...")
        slices_file = kwargs['slices_file']
        constants.N_OF_THREADS=1
        if 'n_of_threads' in kwargs:
            constants.N_OF_THREADS=kwargs['n_of_threads']
        constants.USE_CACHE=False

        if 'use_cache' in kwargs:
            constants.USE_CACHE=kwargs['use_cache']=='true'

        if 'compare_folder' in kwargs:
            compare_folder = kwargs['compare_folder']


        active_genes_file, bg_genes = self.init_params(dataset_file_name, network_file_name, output_folder)
        # print(f'domino_parameters: active_genes_file={active_genes_file}, network_file={network_file_name},slices_file={slices_file}, slice_threshold={slice_threshold},module_threshold={module_threshold}')
        modules = top_sdg_main(dataset_file=dataset_file_name, compare_folder=compare_folder)
        all_bg_genes = [bg_genes for x in modules]
        return modules, all_bg_genes


