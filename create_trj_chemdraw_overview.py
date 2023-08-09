import pandas as pd
import numpy as np
import os
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def SmilesToExactMass(smiles):
    try:

        ps = Chem.SmilesParserParams()
        ps.removeHs = False

        mol = Chem.MolFromSmiles(smiles,ps)

        return rdMolDescriptors.CalcExactMolWt(mol)

    except Exception as e:
        print(f"An error occurred with input {smiles}: {e}")
        return ""

def fill_empty_rows(df_list):
    max_rows = max(df.shape[0] for df in df_list)
    return [df.reindex(pd.RangeIndex(max_rows)) for df in df_list]

def remove_equilibrium_reactions(df):
    df['SMILES'] = df['SMILES'].str.split('.').apply(lambda x: '.'.join([i for i in x if i != 'N#N']))
    df = df.assign(idx = df.index + 1)
    df = df.assign(indx = df.groupby('SMILES').ngroup() + 1)
    df = df.groupby('SMILES').agg({'idx': ['min', 'max'], 'collisions': ['min', 'max']}).reset_index()
    df.columns = ['_'.join(col).strip('_') for col in df.columns.values]
    df = df.sort_values('idx_min').reset_index(drop=True)
    df = df.reset_index(drop=True)
    df = df.assign(id = df.index + 1, flag = False)


    smiles_entry = 1
    while smiles_entry < df.shape[0]:
        df.loc[(df['idx_max'] < df.loc[df['id'] == smiles_entry, 'idx_max'].values[0]) & (df['id'] > smiles_entry), 'flag'] = True
        df = df[df['flag'] == False]
        df = df.reset_index(drop=True)
        df = df.assign(id = df.index + 1)
        smiles_entry += 1

    df = df.groupby(['id', 'SMILES']).apply(lambda group: pd.Series({'collisions': '_'.join(map(str, [group['collisions_min'].values[0], group['collisions_max'].values[0]]))})).reset_index()

    return df

root = 'C:/PostDoc/Ming_time/example_files/csvs_batch'
pattern = '18_protonated_mol_1'

folders = [os.path.join(root, f) for f in os.listdir(root) if pattern in f]
reduce_to_relevant = True

li_allTrj = [None]*len(folders)

pass
for folder_idx, folder in enumerate(folders):

    cid_files = [os.path.join(folder, f) for f in os.listdir(folder) if 'CID' in f]
    MDtrj_files = [os.path.join(folder, f) for f in os.listdir(folder) if 'MDtrj' in f]

    li_singleTrj = [None]*(len(cid_files) + len(MDtrj_files))
    li_idx = 0

    dt = pd.read_csv(MDtrj_files[0], header=None)
    dt['file'] = os.path.basename(MDtrj_files[0])
    dt['collisions'] = 0
    li_singleTrj[li_idx] = dt
    li_idx += 1

    for cid_idx, (cid_file, MDtrj_file) in enumerate(zip(cid_files, MDtrj_files[1:])):

        dt_cid = pd.read_csv(cid_file, header=None)
        dt_md = pd.read_csv(MDtrj_file, header=None)
        dt_cid['file'] = os.path.basename(cid_file)
        dt_md['file'] = os.path.basename(MDtrj_file)
        dt_cid['collisions'] = cid_idx + 1
        dt_md['collisions'] = cid_idx + 1
        li_singleTrj[li_idx*2-1] = dt_cid
        li_singleTrj[li_idx*2] = dt_md
        li_idx += 1

    dt_trj = pd.concat(li_singleTrj, ignore_index=True)

    number = re.search(r'.*TMP\.(\d+)$', folder).group(1)

    if reduce_to_relevant:
        dt_trj = dt_trj.rename(columns={1: 'SMILES'})
        dt_trj = remove_equilibrium_reactions(dt_trj)
        dt_trj.columns = [f'{c}_{number}' for c in ['id', 'SMILES', 'collisions']]
    else:
        dt_trj.columns = [f'{c}_{number}' for c in ['id', 'SMILES', 'file', 'collisions']]

    li_allTrj[folder_idx] = dt_trj

    print(f'{folder_idx + 1}/{len(folders)}')

li_allTrj = fill_empty_rows(li_allTrj)

dt_combined = pd.concat(li_allTrj, axis=1)

dt_combined.to_csv(f'C:/PostDoc/Ming_time/example_files/csvs_summary/summary_' + pattern + '_reduced_python_delete.csv', index=False)


# Iterate over dataframes and modify column names
for i in range(len(li_allTrj)):
    trj_id = li_allTrj[i].columns[0].split('_')[1]
    li_allTrj[i].columns = li_allTrj[i].columns.map(lambda y: y.split('_')[0])
    li_allTrj[i]['trj'] = trj_id

# Combine all dataframes
dt_combinedR = pd.concat(li_allTrj)

# Drop rows where 'id' is NA
dt_combinedR = dt_combinedR.dropna(subset=['id'])

# Add 'smiles_id' column
dt_combinedR['smiles_id'] = dt_combinedR.groupby('SMILES').ngroup()

# Add 'last_smiles' column
def get_smiles_id_with_max_id(group):
    return group.loc[group['id'].idxmax(), 'smiles_id']
dt_combinedR['last_smiles'] = dt_combinedR.groupby('trj').apply(get_smiles_id_with_max_id).reset_index(level=0, drop=True)


# Add 'full_trj' column
# create a dictionary where keys are 'trj' and values are the unique 'smiles_id' in each group
full_trj_dict = dt_combinedR.groupby('trj')['smiles_id'].unique().apply(lambda x: ','.join(map(str, x))).to_dict()

# map 'trj' column to 'full_trj' values using the dictionary
dt_combinedR['full_trj'] = dt_combinedR['trj'].map(full_trj_dict)


dt_combinedR['mz'] = dt_combinedR['SMILES'].apply(SmilesToExactMass)




# Generate 'dt_full_trj_summary'
# Group by 'trj_id' and 'trj', and calculate the size of each group

grouped = dt_combinedR.groupby('full_trj')['trj'].nunique().reset_index(name='n')


# Merge back with the original DataFrame to get 'SMILES' and 'mz' columns
# Get the first 'trj' value of each 'full_trj' group
first_trj_per_group = dt_combinedR.groupby('full_trj')['trj'].transform('first')

# Create a boolean mask that is True for rows where 'trj' is the first 'trj' of its 'full_trj' group
mask = dt_combinedR['trj'] == first_trj_per_group

# Use the mask to filter the DataFrame
dt_combinedR_tmp = dt_combinedR[mask]


dt_full_trj_summary = pd.merge(grouped, dt_combinedR_tmp[['id','full_trj', 'SMILES', 'mz']], on='full_trj', how='left')



dt_combinedR.to_csv(f'C:/PostDoc/Ming_time/example_files/csvs_summary/summary_' + pattern + '_dt_combinedR.csv', index=False)
dt_full_trj_summary.to_csv(f'C:/PostDoc/Ming_time/example_files/csvs_summary/summary_' + pattern + '_dt_full_trj_summary.csv', index=False)