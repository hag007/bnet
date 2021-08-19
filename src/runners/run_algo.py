import sys
sys.path.insert(0, '../')
import os
from src import constants
from src.runners.bnet_runner import BnetRunner
from src.runners.domino_runner import DominoRunner
from src.runners.bnet_static_string_runner import BnetStaticStringRunner
from src.runners.bnet_static_corr_runner import BnetStaticCorrRunner
from src.runners.netbox_runner import NetboxRunner
from src.runners.jactivemodules_greedy_runner import jAMGreedyRunner
from src.runners.jactivemodules_sa_runner import jAMSARunner
from src.runners.bionet_runner import BionetRunner
from src.runners.keypathwayminer_ines_greedy_runner import KPMRunner
from src.runners.hotnet2_runner import Hotnet2Runner

ALGO_BY_NAMES = {"netbox": NetboxRunner(), "jactivemodules_greedy": jAMGreedyRunner(),
                 "jactivemodules_greedy_string": jAMGreedyRunner(), "jactivemodules_sa": jAMSARunner(),
                 "bionet": BionetRunner(), "bionet_string": BionetRunner(), 'keypathwayminer_INES_GREEDY': KPMRunner(),
                 'hotnet2': Hotnet2Runner(), 'BNET_STATIC_STRING': BnetStaticStringRunner(), 'BNET_STATIC_CORR': BnetStaticCorrRunner(),
                 'BNET_dynamic_Cosine': BnetRunner('Cosine'), 'BNET_dynamic_Dice': BnetRunner('Dice'), 'BNET_dynamic_Jaccard': BnetRunner('Jaccard'),
                 'BNET_dynamic_Lin': BnetRunner('Lin'), 'BNET_dynamic_Jiang-Conrath': BnetRunner('Jiang-Conrath'),
                 'BNET_dynamic_SimGIC': BnetRunner('SimGIC'), 'BNET_dynamic_SimRel': BnetRunner('SimRel'), 'BNET_dynamic_SimUI': BnetRunner('SimUI'),
                 'top_SDG_genes': TopSDGGenesRunner('top_SDG_genes'), 'top_SDG_pcs': TopSDGPcsRunner('top_SDG_pcs'), 'DOMINO': BnetRunner('DOMINO')} #




def add_algo_runner(k,v):
    ALGO_BY_NAMES[k]=v

def create_ds_folders(dataset_name):
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "data")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "cache")))
    os.makedirs(os.path.join(os.path.join(constants.DATASETS_DIR, dataset_name, "output")))


def run_algo(dataset_name, algo, network_file_name, go_folder, output_folder, **kwargs):
    ALGO_BY_NAMES[algo].main(dataset_name, network_file_name, go_folder, output_folder, **kwargs)
#
#
