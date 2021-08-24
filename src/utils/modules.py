import os
import pandas as pd

def get_num_modules(compare_folder):
    df = pd.read_csv(os.path.join(compare_folder, "modules/modules_summary.tsv"))
    return df.shape[0]