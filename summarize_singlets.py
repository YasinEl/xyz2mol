import pandas as pd
import numpy as np
import os
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import argparse
import json
import glob


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

        li_singleTrj = []

        dt = pd.read_csv(MDtrj_files[0], header=None)
        dt['file'] = os.path.basename(MDtrj_files[0])
        dt['collisions'] = 0
        li_singleTrj.append(dt)

        li_singleTrj = []
        cid_files = sorted(cid_files, key=lambda x: int(x.split("CID")[1].split('.')[0]))
        MDtrj_files = sorted(MDtrj_files[1:], key=lambda x: (int(x.split(".")[-3]), int(x.split(".")[-2])))

        for cid_file in cid_files:
            dt_cid = pd.read_csv(cid_file, header=None)
            dt_cid['file'] = os.path.basename(cid_file)
            cid_idx = int(cid_file.split("CID")[1].split('.')[0])
            dt_cid['collisions'] = cid_idx
            li_singleTrj.append(dt_cid)

            # Filter and sort the corresponding MDtrj files for the current CID file
            corresponding_MDtrj_files = sorted([f for f in MDtrj_files if f.split(".")[-3] == str(cid_idx)],
                                               key=lambda x: int(x.split(".")[-2]))

            for MDtrj_file in corresponding_MDtrj_files:
                dt_md = pd.read_csv(MDtrj_file, header=None)
                dt_md['file'] = os.path.basename(MDtrj_file)
                dt_md['collisions'] = cid_idx
                li_singleTrj.append(dt_md)

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

def list_relevant_csv_files(directory, pattern):

    relevant_files = []

    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith('.csv') and (pattern + '__' in filename): #change for numbers
                relevant_files.append(os.path.join(dirpath, filename))

    relevant_files = [f for f in relevant_files if re.match(r'^[\d.]+\.\w+$', f.split('__')[-1])]

    return relevant_files



def summarize_singlets(outfile_json, structure_csvs):
    filename = os.path.basename(outfile_json).split(".json")[0]
    parts = filename.split("__")[0]

    with open(outfile_json, 'r') as file:
        data = json.load(file)

    csv_files = sorted(list_relevant_csv_files(structure_csvs, parts))

    # Prepare a dataframe to hold all combined data
    df_list = []

    # Loop through each csv file, read its data and add it to the main dataframe
    for file in csv_files:
        base = os.path.basename(file).rsplit('.', 1)[0]
        collision, run_num, fragment_num = base.split(".tmol_TMP.")[1].split("__")[1].split(".")

        temp_df = pd.read_csv(file, header=None, names=['ID', 'SMILES'])
        temp_df['Source CSV'] = os.path.basename(file)
        temp_df['Collision'] = collision
        temp_df['Run-Number'] = run_num
        temp_df['Fragment Number'] = fragment_num

        df_list.append(temp_df)

    # Combine all dataframes row-wise
    main_df = pd.concat(df_list, ignore_index=True)

    # Add information from json based on the csv filename
    for idx, row in main_df.iterrows():
        apply_offset = 0
        add_trj_offset = 0
        identifier = ""
        if row['Collision'] == "0":
            identifier = "heating_trajectory" + row['Run-Number']
        elif row['Run-Number'] == "1":
            identifier = "CID" + row['Collision']
        elif int(row['Run-Number']) > 1:
            identifier = "MD" + row['Collision']
            if int(row['Run-Number']) > 2:
                add_trj_offset = int(row['Run-Number']) - 1
        else:
            continue

        # Check if the identifier exists in the data
        if identifier not in data:
            print(f"Warning: Identifier '{identifier}' not found in JSON data.")
            continue

        # Check if the ID matches
        if int(row['Fragment Number']) in [charge['id'] for charge in data[identifier]['charges_per_fragment']]:
            if add_trj_offset > 0:
                count_id1 = 0
                for idx2, item in enumerate(data[identifier]['charges_per_fragment']):
                    if item['id'] == 1:
                        count_id1 += 1
                    if count_id1 == add_trj_offset:
                        break
                    apply_offset += 1

            main_df.at[idx, 'diss time (ps)'] = data[identifier]['fragments']['diss time (ps)'][int(row['Fragment Number']) - 1 + apply_offset]
            main_df.at[idx, 'formula'] = data[identifier]['fragments']['formula'][int(row['Fragment Number']) - 1  + apply_offset]
            main_df.at[idx, 'mass'] = data[identifier]['fragments']['mass'][int(row['Fragment Number']) - 1  + apply_offset]
            main_df.at[idx, 'charge'] = data[identifier]['charges_per_fragment'][int(row['Fragment Number']) - 1  + apply_offset]['charge']
            # Add other required fields in a similar manner

            main_df.to_csv(parts + '__SingletsSummary.csv')

    return main_df



if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Summarize single fragments')
    parser.add_argument('--input_json', '-ij',type=str, required=True, help='Path to json file')
    parser.add_argument('--input_csvs', '-ic', type=str, required=True, help='Path to csv files')
    args = parser.parse_args()

    df = summarize_singlets(args.input_json, args.input_csvs)
