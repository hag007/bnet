import pandas as pd
import os
file = "/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split.tsv"

file_df = pd.read_csv(file, sep='\t')
train_files = ["bnet_brca_M","bnet_brca_N","bnet_brca_stage","bnet_brca_T","bnet_luad_retrospective","bnet_luad_source","bnet_luad_vital","bnet_ov_source_site",
               "bnet_ov_vital","bnet_skcm_location","bnet_skcm_retrospective","bnet_skcm_status","bnet_ucec_hypertension","bnet_ucec_retrospective","bnet_ucec_type","bnet_ucec_vital"]
test_files = []


train_df = pd.DataFrame(columns=file_df.columns)
test_df = pd.DataFrame(columns=file_df.columns)

for i in range(len(file_df)):
    name = file_df.iloc[i]['path to genes']
    relevant = name.split('/')[-1].split('.')[0]
    if relevant in train_files:
        train_df = train_df.append(file_df.iloc[i], ignore_index=True)
    else:
        test_df = test_df.append(file_df.iloc[i], ignore_index=True)

train_df.to_csv("/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split_real_train.tsv", sep='\t', index=False, header=False)
os.chmod("/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split_real_train.tsv", 0o777)
test_df.to_csv("/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split_real_test.tsv", sep='\t', index=False, header=False)
os.chmod("/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split_real_test.tsv", 0o777)

