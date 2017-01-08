import os
import pandas as pd

RAW_PATH = './reference_data'

with open('references.tsv', 'w') as f:
    for filename in os.listdir(RAW_PATH):
        file_obj = open(os.path.join(RAW_PATH, filename), 'rb')
        raw_table = pd.read_csv(file_obj, sep='\t')

        f.write(filename + '\t')
        f.write(','.join(set(raw_table['Ref Loc'])))
        f.write('\n')
