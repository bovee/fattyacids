# run this over the references to check that all their names are
# within NCBI
import os
import pandas as pd

RAW_PATH = './reference_data'

species_ids = set()
genus_ids = set()
parent_id = {}
with open('../taxonomy/taxdmp/nodes.dmp', 'r') as nodes:
    for n in nodes:
        line = n.split('\t|\t')
        parent_id[line[0]] = line[1]
        if line[2] == 'species':
            species_ids.add(line[0])
        elif line[2] == 'genus':
            genus_ids.add(line[0])


name_to_id = {}
with open('../taxonomy/taxdmp/names.dmp', 'r') as names:
    for l in names:
        tax_id, name, _ = l.split('\t|\t', 2)
        if tax_id in genus_ids or tax_id in species_ids:
            name_to_id[name] = tax_id


def species_to_taxid(filename, name):
    if len(name.split()) > 2:
        if name.split()[2] in {'forma', 'subspecies', 'serovar'}:
            name = ' '.join(name.split()[:2])
        elif name.split()[1] == 'cf':
            # maybe we should drop these entirely?
            name = name.split()[0]
    if name not in name_to_id:
        if name.startswith('“') or name == '':
            # print(filename, name)
            pass
        else:
            print(filename, name)
        return ''
    else:
        return name_to_id[name]


for filename in os.listdir(RAW_PATH):
    file_obj = open(os.path.join(RAW_PATH, filename), 'rb')
    raw_table = pd.read_csv(file_obj, sep='\t')

    # del raw_table['Ref Loc']
    for s in raw_table['Species'].fillna(''):
        species_to_taxid(filename, s)
    # raw_table.insert(0, 'Tax ID', tax_col)
    # raw_table.to_csv(os.path.join(RAW_PATH, filename), sep='\t', index=False)
