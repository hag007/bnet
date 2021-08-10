import os
import pandas as pd
import json

df = pd.DataFrame(columns=["path to genes", "path to TCGA file", "name of variable", "value1", "value2"])


directory_str = "/specific/netapp5/gaga/hagailevi/omics/input"
directory = os.fsencode(directory_str)
index = 0
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.startswith("bnet_"):
        genes_file = directory_str+"/"+filename
        project = filename.split("_")[1].upper()
        tcga_file = "/specific/netapp5/gaga/hagailevi/omics/GDC-TCGA/"+project+"/tcga_data/TCGA-"+project+".GDC_phenotype.tsv"
        json_file = "/specific/netapp5/gaga/hagailevi/melanoma/groups/"+filename.split(".")[0]+".json"
        f = open(json_file, )
        data = json.load(f)
        key = list(data[0].keys())[0]
        val1 = data[0][key]['value']
        val2 = data[1][key]['value']
        lst = [genes_file, tcga_file, key, val1, val2]
        df.loc[index] = lst
        index+=1
    else:
        continue
df.to_csv("/specific/netapp5/gaga/hagailevi/omics/input/phenotypes_split.tsv", sep="\t", index=False)
print(0)