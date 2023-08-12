import pandas as pd
import numpy as np
import os
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import argparse
import json
import glob
from molmass import Formula


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
            main_df.at[idx, 'formula'] = data[identifier]['fragments']['formula'][int(row['Fragment Number']) - 1 + apply_offset]
            main_df.at[idx, 'mass'] = data[identifier]['fragments']['mass'][int(row['Fragment Number']) - 1 + apply_offset]
            main_df.at[idx, 'mass_from_formula'] = Formula(data[identifier]['fragments']['formula'][int(row['Fragment Number']) - 1 + apply_offset] + "+").isotope.mass
            main_df.at[idx, 'charge'] = data[identifier]['charges_per_fragment'][int(row['Fragment Number']) - 1 + apply_offset]['charge']
            # Add other required fields in a similar manner

    return main_df, parts


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Get fragment edges')
    parser.add_argument('--singlet_csvs', '-ij',type=str, required=True, help='Path to json file')
    parser.add_argument('--input_csvs', '-ic', type=str, required=True, help='Path to csv files')
    args = parser.parse_args()

    df, trj_name = summarize_singlets(args.input_json, args.input_csvs)

    remove_equilibrium_reactions(df)

    df_edges = df[['SMILES']]

    df_edges['target'] = df_edges['SMILES'].shift(-1)
    df_edges['source'] = df_edges['SMILES']

    df_edges = df_edges.drop_duplicates(subset='target', keep='last')

    df.to_csv(trj_name + '__SingletsSummary.csv')
    df_edges.to_csv(trj_name + '__SingletsEdges.csv')

