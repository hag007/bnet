
import pandas as pd
import numpy as np

cancer_types=["BRCA","OV","SKCM","LUAD","UCEC"]
network_files=["dip","huri","string"]

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
        df_network.loc[:,"score"]=df_network.apply(lambda a: f"{{'weight': {mean if not a[0] in df_corr.index or not a[2] in df_corr.index else df_corr.loc[a[0], a[2]]} }}",axis=1)
        df_network.iloc[:,[0,2,3]].to_csv(f'/specific/netapp5/gaga/hagailevi/omics/input/{network_file}_{cancer_type}.txt', sep='\t', index=False)
