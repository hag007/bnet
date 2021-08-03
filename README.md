# Bnet


## Outline


- [Set your environment](#set-your-environment)
- [Run Bnet](#run-bnet)
- [Main output files](#main-output-files)
- [Bnet container](#bnet-container)

## Set your environment

Download the sources and install according to the following instruction:

Clone the repo from github:
```
git clone https://github.com/hag007/bnet.git
cd bnet
```

Bnet is written in Python 3.6. We recommend using a virtual environment. in Linux:
```
python3 -m venv bnet-env
source bnet-env/bin/activate
```

To install Bnet dependencies type:
```
pip install -r  config/dependencies.txt
```

## Run Bnet

Bnet consists of several steps. For a specific set of input parameters, these steps should be carried sequentially.  
Each parameter can be specified as command line parameter (For example `python script.py --param1 value1 --param2 value2`). values of parameters which are not specified in the command line are taken from `config/conf.json`.      

  
1. `generate_solution.py`: Run the bnet_sa algorithm.  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: the algorithm to execute.  
`--permuted_solutions_folder`: folder where permuted solutions reside.  
`--true_solutions_folder`: folder where true solutions reside.  
`--go_folder`: folder where GO files are located.  
`--network_file`: file of the biological network of the analysis.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm. 

2. `calc_pcs.py`: For each module reported in the bnet_sa solution, the first PC is extracted. These PCs will server later as features to train classifiers (see step #3) .  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: the algoritm to execute.  
`--network_file`: file of the biological network of the analysis.  
`--go_folder`: folder where GO files are located.  
`--true_solutions_folder`: folder where true solutions reside.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm. 

3. `calc_prediction.py`: uses the PCs generated in step #2 as features to train the classifiers RF and SVM. In addition, for each classifier, and check the following metrics: F1, AUPR, AUROC.  
parameters:  
`--dataset_file`: path to dataset file.  
`--algo`: AMI algorithm.  
`--network_file`: file of the biological network of the analysis.  
`--go_folder`: folder where GO files are located.  
`--true_solutions_folder`: folder where true solutions reside.  
`--additional_args`: additional arguments that are relevant to a particular AMI algorithm. 

Parameters default values are defined at `config/conf.json`  

## Main output files

TBD

## Bnet container
Bnet is also available as ready-to-use tool in a container.
Using TAU-VPN and gdocker, type `gdocker gaga-import --image_name=bnet`

