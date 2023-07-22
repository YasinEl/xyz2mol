import pandas as pd
import numpy as np
import os
import re

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
    df = df.assign(id = df.index + 1, flag = False)

    smiles_entry = 1
    while smiles_entry < df.shape[0]:
        df.loc[(df['idx_max'] < df.loc[df['id'] == smiles_entry, 'idx_max'].values[0]) & (df['id'] > smiles_entry), 'flag'] = True
        df = df[df['flag'] == False]
        df = df.assign(id = df.index + 1)
        smiles_entry += 1

    df = df.groupby(['id', 'SMILES']).apply(lambda group: pd.Series({'collisions': '_'.join(map(str, [group['collisions_min'].values[0], group['collisions_max'].values[0]]))})).reset_index()

    return df

root = 'C:/PostDoc/Ming_time/example_files/csvs'

folders = [os.path.join(root, f) for f in os.listdir(root) if '2_protonated_mol_3' in f]
reduce_to_relevant = True

li_allTrj = [None]*len(folders)

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
        dt = pd.read_csv(cid_file if 'CID' in cid_file else MDtrj_file, header=None)
        dt['file'] = os.path.basename(cid_file if 'CID' in cid_file else MDtrj_file)
        dt['collisions'] = cid_idx
        li_singleTrj[li_idx] = dt
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

dt_combined.to_csv(f'C:/PostDoc/Ming_time/example_files/csvs_summary/summary_2_protonated_mol_3_reduced_python.csv', index=False)
