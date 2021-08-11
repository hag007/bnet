is_greedy=True
algo_dir=/specific/netapp5/gaga/hagailevi/evaluation/bnetworks_alg/jactivemodules
num_of_modules=50
overlap_threshold=0
network_file_name=/specific/netapp5/gaga/hagailevi/EMP/data/emp_test/networks/dip.sif
score_file_name=/specific/netapp5/gaga/hagailevi/EMP/data/emp_test/original_datasets/tnfa.tsv
output_file=/specific/netapp5/gaga/hagailevi/EMP/data/emp_test/true_solutions/tnfa_jactivemodules_greedy/jactivemodules_greedy_results.txt


echo /specific/netapp5/gaga/hagailevi/EMP/data/emp_test/original_datasets/tnfa.tsv
java -jar $algo_dir/jactivemodules.jar \
                                 $network_file_name \
                                 $score_file_name \
                                 $is_greedy \
                                 $num_of_modules \
                                 $overlap_threshold \
                                 $output_file