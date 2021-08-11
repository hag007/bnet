
import pandas as pd
import numpy as np
import os

cancer_types=["BRCA","OV","SKCM","LUAD","UCEC"]
network_files=["dip","huri","string"]

string_file = "/specific/netapp5/gaga/hagailevi/omics/input/string_genes_all.txt"
df_string = pd.read_csv(string_file, sep='\t', index_col=0)
df_string['index1'] = df_string.index

string_dict = {}
for i in range(len(df_string)):
# for i in range(100):
    key = frozenset([df_string.iloc[i][0], df_string.iloc[i][2]])
    string_dict[key]=df_string.iloc[i][1]

keys = string_dict.keys()
for network_file in network_files:
    string_network_df = pd.DataFrame(columns=["p1", "p2", "weight"])
    not_sound_pp = pd.DataFrame()
    df_network = pd.read_csv(f'/specific/netapp5/gaga/hagailevi/omics/input/{network_file}.sif', sep='\t')
    index = 0
    for i in range(len(df_network)):
        pp = frozenset([df_network.iloc[i][0], df_network.iloc[i][2]])
        if pp in keys:
            w = "{\'score\':"+str(string_dict[pp])+"}"
            string_network_df.loc[index] = [df_network.iloc[i][0], df_network.iloc[i][2], w]
            index += 1
        else:
            content=list(pp)
            if len(content)==1:
                content.append(content[0])
            not_sound_pp = not_sound_pp.append({"p0":content[0], "p1":content[1]}, ignore_index=True)
    file_name = "/specific/netapp5/gaga/hagailevi/omics/input/string_"+network_file+".tsv"
    no_found_file_name = "/specific/netapp5/gaga/hagailevi/omics/input/not_found_string_"+network_file+".tsv"
    string_network_df.to_csv(file_name, sep="\t", index=False, header=False, )
    not_sound_pp.to_csv(no_found_file_name, sep="\t", index=False, header=False)
    os.chmod(file_name, 0o777)
    os.chmod(no_found_file_name, 0o777)
    print("done!")

#
for cancer_type in cancer_types:
    for network_file in network_files:
        df=pd.read_csv(f'/specific/netapp5/gaga/hagailevi/omics/GDC-TCGA/{cancer_type}/tcga_data/TCGA-{cancer_type}.htseq_fpkm.tsv', sep='\t', index_col=0)
        df_index=pd.read_csv('/specific/netapp5/gaga/hagailevi/omics/list/protein_coding.txt', sep='\t', index_col=0).index
        df.index=[a.split(".")[0] for a in df.index]
        df_clean=df.reindex(df_index).dropna(axis=0)
        res=np.corrcoef(df_clean)
        df_corr=pd.DataFrame(index=df_clean.index, columns=df_clean.index,data=res)

        df_network=pd.read_csv(f'/specific/netapp5/gaga/hagailevi/omics/input/{network_file}.sif', sep='\t')
        mean = df_network.apply(lambda a: None if not a[0] in df_corr.index or not a[2] in df_corr.index else df_corr.loc[a[0], a[2]],axis=1).mean()
        df_network.loc[:,"score"]=df_network.apply(lambda a: f"{{'score': {mean if not a[0] in df_corr.index or not a[2] in df_corr.index or np.isnan(df_corr.loc[a[0], a[2]]) else df_corr.loc[a[0], a[2]]} }}",axis=1)
        df_network.iloc[:,[0,2,3]].to_csv(f'/specific/netapp5/gaga/hagailevi/omics/input/{network_file}_{cancer_type}.txt', sep='\t', index=False, header=False)
        os.chmod(f'/specific/netapp5/gaga/hagailevi/omics/input/{network_file}_{cancer_type}.txt', 0o777)
