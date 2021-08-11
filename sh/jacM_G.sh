is_greedy=True
algo_dir=/specific/elkon/hagailevi/DOMINO_FILES/evaluation/bnetworks_alg/jactivemodules
num_of_modules=5
overlap_threshold=0
network_file_name=network2.sif
score_file_name=test2.tsv
output_file="/specific/netapp5/gaga/hagailevi/EMP/src/sh/out.txt"


echo $score_file_name
cmd="java -jar $algo_dir/jam.jar \
                                 $network_file_name \
                                 $score_file_name \
                                 $is_greedy \
                                 $num_of_modules \
                                 $overlap_threshold \
                                 $output_file"


echo $cmd

eval $cmd

