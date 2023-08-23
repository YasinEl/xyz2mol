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
            if filename.endswith('.csv') and (pattern + '__' in filename): 
                relevant_files.append(os.path.join(dirpath, filename))

    relevant_files = [f for f in relevant_files if re.match(r'^[\d.]+\.\w+$', f.split('__')[-1]) or f.split('__')[-1] == 'start.csv']

    def custom_sort(item):
        if "start.csv" in item:
            return ("", item)
        return (item,)

    return sorted(relevant_files, key = custom_sort) 



def summarize_singlets(outfile_json, structure_csvs):
    filename = os.path.basename(outfile_json).split(".json")[0]
    parts = filename.split("__")[0]

    with open(outfile_json, 'r') as file:
        data = json.load(file)

    if data['status'] == 'failed':
        return 'failed', 'failed'

    csv_files = list_relevant_csv_files(structure_csvs, parts)


    if len(csv_files) == 0:
        return 'failed', 'failed'
    
    # Prepare a dataframe to hold all combined data
    df_list = []

    # Loop through each csv file, read its data and add it to the main dataframe
    for file in csv_files:
        base = os.path.basename(file).rsplit('.', 1)[0]
        if not base.endswith('start'):
            collision, run_num, fragment_num = base.split(".tmol_TMP.")[1].split("__")[1].split(".")
        else:
            collision = 0
            run_num = 0
            fragment_num = 0


        temp_df = pd.read_csv(file, header=None, names=['ID', 'SMILES'])
        temp_df['Source CSV'] = os.path.basename(file)
        temp_df['Collision'] = int(collision)
        temp_df['Run-Number'] = int(run_num)
        temp_df['Fragment Number'] = int(fragment_num)

        df_list.append(temp_df)

    # Combine all dataframes row-wise
    main_df = pd.concat(df_list, ignore_index=True)

    # Add information from json based on the csv filename
    for idx, row in main_df.iterrows():
        identifier = ""
        if row['Collision'] == 0:
            identifier = "heating_trajectory" + str(row['Run-Number'])
        elif row['Run-Number'] == 1:
            identifier = "CID" + str(row['Collision'])
        elif row['Run-Number'] > 1:
            identifier = "MD" + str(row['Collision']) + '_' + str(row['Run-Number'] - 1)

        else:
            continue

        # Check if the identifier exists in the data
        if identifier not in data:
            print(f"Warning: Identifier '{identifier}' not found in JSON data.")
            continue

        # Check if the ID matches
        if row['Fragment Number'] in [charge['id'] for charge in data[identifier]['charges_per_fragment']]:

            main_df.at[idx, 'diss time (ps)'] = data[identifier]['fragments']['diss time (ps)'][row['Fragment Number'] - 1]
            main_df.at[idx, 'formula'] = data[identifier]['fragments']['formula'][row['Fragment Number'] - 1]
            main_df.at[idx, 'mass'] = data[identifier]['fragments']['mass'][row['Fragment Number'] - 1]
            main_df.at[idx, 'mass_from_formula'] = Formula(data[identifier]['fragments']['formula'][row['Fragment Number'] - 1] + "+").isotope.mass
            main_df.at[idx, 'charge'] = data[identifier]['charges_per_fragment'][row['Fragment Number'] - 1]['charge']
            main_df.at[idx, 'nominal_charge'] = 1 if data[identifier]['charges_per_fragment'][row['Fragment Number'] - 1]['used'] == True else 0

            # Add other required fields in a similar manner

    #main_df['nominal_charge'] = main_df.groupby(['Collision', 'Run-Number'])['charge'].transform(lambda x: x == x.max()).astype(int)
    main_df.loc[main_df['Source CSV'].str.contains('start.csv'), 'nominal_charge'] = 1

    main_df['Collision'] = pd.to_numeric(main_df['Collision'], errors='coerce')
    main_df['Run-Number'] = pd.to_numeric(main_df['Run-Number'], errors='coerce')

    main_df = main_df.sort_values(by=['Collision', 'Run-Number'])

    return main_df, parts


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Get fragment edges')
    parser.add_argument('--input_json', '-ij',type=str, required=True, help='Path to json file')
    parser.add_argument('--input_csvs', '-ic', type=str, required=True, help='Path to csv files')
    args = parser.parse_args()


    df, trj_name = summarize_singlets(args.input_json, args.input_csvs)

    if trj_name != 'failed':

        df.to_csv(trj_name + '__SingletsSummary.csv', index=False)

        #get edges table
        df_edges = pd.DataFrame(columns=['source', 'target'])

        prec = None
        prod = None
        next_prec = None

        for index, row in df.iterrows():

            if prec is None:
                next_prec = row['SMILES']

            if prec is None or (row['Collision'] > current_col or row['Run-Number'] > current_run):
                prec = next_prec

            current_smiles = row['SMILES']
            if row['SMILES'] != prec:
                prod = row['SMILES']

            if len(df_edges[(df_edges['source'] == prec) & (df_edges['target'] == prod)]) == 0 and not prod is None and prec != prod:
                new_row = pd.DataFrame({
                    'source': [prec],
                    'target': [prod]
                })
                df_edges = pd.concat([df_edges, new_row], ignore_index=True)

            if prec is None or row['nominal_charge'] == 1:
                current_col = row['Collision']
                current_run = row['Run-Number']
                next_prec = row['SMILES']


        if len(df_edges) > 0:
            df_edges.to_csv(trj_name + '__SingletsEdges.csv', index=False)

