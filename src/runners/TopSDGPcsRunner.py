import sys
sys.path.insert(0, '../')
sys.path.insert(0, '../..')
import os
import pandas as pd

from src import constants
from src.implementations.PCs_as_genes import main as PCs_as_genes_main
from src.utils.ensembl2entrez import ensembl2entrez_convertor
from src.utils.network import get_network_genes
from src.utils.go_similarity import init_go_metadata

from src.runners.abstract_runner import AbstractRunner


class PCsAsNoGenesRunner(AbstractRunner):
    def __init__(self):
        super().__init__(f"PCs_as_genes")

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
        print("run bnet_dynamic runner...")
        slices_file = kwargs['slices_file']
        constants.N_OF_THREADS=1
        if 'n_of_threads' in kwargs:
            constants.N_OF_THREADS=kwargs['n_of_threads']
        constants.USE_CACHE=False
        if 'use_cache' in kwargs:
            constants.USE_CACHE=kwargs['use_cache']=='true'
        slice_threshold = 0.3
        if 'slice_threshold' in kwargs:
            slice_threshold = kwargs['slice_threshold']
        module_threshold = 0.05
        if 'module_threshold' in kwargs:
            module_threshold = kwargs['module_threshold']
        algo = "BNET_STATIC_STRING"
        if 'algo' in kwargs:
            algo = kwargs['algo']
        true_solutions_folder = "/home/gaga/hagailevi/omics/output/true_solutions"
        if 'true_solutions_folder' in kwargs:
            true_solutions_folder = kwargs['true_solutions_folder']
        ts = 100
        if 'ts' in kwargs:
            ts = kwargs['ts']
        min_temp = 10
        if 'min_temp' in kwargs:
            min_temp = kwargs['min_temp']
        temp_factor = 40.0
        if 'temp_factor' in kwargs:
            temp_factor = kwargs['temp_factor']
        qval_norm = 1.3
        if 'qval_norm' in kwargs:
            qval_norm = kwargs['qval_norm']
        min_n_genes = 4
        if 'min_n_genes' in kwargs:
            min_n_genes = kwargs['min_n_genes']
        sim_factor = 2.5
        if 'sim_factor' in kwargs:
            sim_factor = kwargs['sim_factor']
        activity_baseline = 0
        if 'activity_baseline' in kwargs:
            activity_baseline = kwargs['activity_baseline']

        active_genes_file, bg_genes = self.init_params(dataset_file_name, network_file_name, output_folder)
        modules = PCs_as_genes_main(dataset_file=dataset_file_name, network_file=network_file_name,
                            slice_threshold=slice_threshold, module_threshold=module_threshold, algo=algo,
                               true_solutions_folder=true_solutions_folder, ts=ts, min_temp=min_temp, temp_factor=temp_factor,
                             sim_factor=sim_factor, activity_baseline=activity_baseline)
        all_bg_genes = [bg_genes for x in modules]
        return modules, all_bg_genes


