import pandas as pd
import numpy as np
import os
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import argparse

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


def list_relevant_csv_files(directory):
    """
    List all .csv files with 'CID' or 'MDtrj' in their names
    from the given directory and its subdirectories.
    """
    relevant_files = []

    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.csv') and ('CID' in filename or 'MDtrj' in filename):
                relevant_files.append(os.path.join(dirpath, filename))

    return relevant_files



def summarize_trajectories(input_directory):
    csv_files = list_relevant_csv_files(input_directory)
    unique_trajectories = set()

    for csv_file in csv_files:
        filename = os.path.basename(csv_file).split(".csv")[0]
        parts = filename.split("__")
        if len(parts) == 2:
            unique_trajectories.add(parts[0])

    li_allTrj = [None] * len(unique_trajectories)
    all_files = os.listdir(input_directory)

    for trajectory_idx, trajectory_name in enumerate(unique_trajectories):
        cid_files = [os.path.join(input_directory, f) for f in all_files if trajectory_name in f and 'CID' in f]
        MDtrj_files = [os.path.join(input_directory, f) for f in all_files if trajectory_name in f and 'MDtrj' in f]


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

        number = re.search(r'.*TMP\.(\d+)$', trajectory_name).group(1)

        #reduce to relevant
        dt_trj = dt_trj.rename(columns={1: 'SMILES'})
        dt_trj = remove_equilibrium_reactions(dt_trj)
        dt_trj.columns = [f'{c}_{number}' for c in ['id', 'SMILES', 'collisions']]

        li_allTrj[trajectory_idx] = dt_trj

    li_allTrj = fill_empty_rows(li_allTrj)

    dt_combined = pd.concat(li_allTrj, axis=1)

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

    group_reprs = dt_combinedR.groupby('trj').apply(
        lambda group: ''.join([str(x) + '_' + str(y) for x, y in zip(group['id'], group['smiles_id'])])).reset_index(
        name='group_repr')

    # Calculate the count for each group_repr
    group_reprs_counts = group_reprs['group_repr'].value_counts().reset_index()

    group_reprs_counts.columns = ['group_repr', 'count']

    # Merge the count back to the group_reprs dataframe
    group_reprs = pd.merge(group_reprs, group_reprs_counts, on='group_repr')

    # Merge the group_repr back to the original dataframe on 'trj'
    dt_combinedR = pd.merge(dt_combinedR, group_reprs, on='trj')

    def concatenate_unique(x):
        # Convert the series to set (for unique values), then back to list, then join as string
        return ', '.join(map(str, sorted(set(x))))

    # Group by 'group_repr' and aggregate the 'trj' column using the above function
    df_unique = dt_combinedR.groupby('group_repr')['trj'].agg(concatenate_unique).reset_index()

    # Merge the resulting DataFrame with the original one on 'group_repr' to add the new column
    dt_combinedR = dt_combinedR.merge(df_unique, on='group_repr', how='left').rename(columns={'trj_x': 'trj', 'trj_y': 'unique_trj_values'})


    # Filter to keep only one representative group for each set of identical groups
    dt_combinedR = dt_combinedR.drop_duplicates(subset=['smiles_id', 'group_repr']).drop(columns=['group_repr', 'collisions', 'trj'])
    dt_combinedR['mz'] = dt_combinedR['SMILES'].apply(SmilesToExactMass)


    dt_combined.to_csv(os.path.join("collectedTrajectories.csv"), index=False)
    dt_combinedR.to_csv(os.path.join("overall_rows.csv"), index=False)

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Summarize trajectories')
    parser.add_argument('--input', '-i',type=str, required=True, help='Path to csv files')
    args = parser.parse_args()

    summarize_trajectories(args.input)

    #summarize_trajectories("C:/Users/elabi/Downloads/csvs", "C:/PostDoc/Ming_time/example_files/test")
