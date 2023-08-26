import pandas as pd
import os
import glob
import argparse

def main(input_path):
    # 1st CSV processing
    edge_files = glob.glob(os.path.join(input_path, '*__SingletsEdges.csv'))
    dfs = [pd.read_csv(file) for file in edge_files]
    dt = pd.concat(dfs, axis=0, ignore_index=True)
    dt = dt.groupby(['source', 'target']).size().reset_index(name='n')
    dt.dropna(subset=['source'], inplace=True)
    dt.to_csv("Edge_table.csv", index=False)

    # 2nd CSV processing
    summary_files = glob.glob(os.path.join(input_path, '*__SingletsSummary.csv'))
    dfs = [pd.read_csv(file) for file in summary_files]
    dt = pd.concat(dfs, axis=0, ignore_index=True)
    dt['trj'] = dt['Source CSV'].str.split("__").str[0]
    dt['ID'] = dt.groupby('trj').cumcount() + 1

    dt['nominal_charge_temp'] = dt.groupby(['trj', 'Collision', 'Run-Number'])['nominal_charge'].transform(max)
    max_charge = dt[dt['nominal_charge_temp'] == 0]['charge'].max()
    dt.loc[dt['nominal_charge_temp'] == 0, 'nominal_charge'] = max_charge

    dt = dt.drop(columns='nominal_charge_temp')

    dt['id_max_charged'] = dt[dt['nominal_charge'] == 1].groupby('trj')['ID'].transform(max)
    dt['trj_end'] = dt['id_max_charged'] == dt['ID']

    agg_funcs = {
        'mass_from_formula': lambda x: ','.join(map(str, x[(x != '') & pd.notna(x)].unique())),
        'nominal_charge': 'mean',
        'Collision': ['min', 'max', 'mean'],
        'trj_end': 'any',
        'formula': lambda x: ','.join(map(str, x[(x != '') & pd.notna(x)].unique()))
    }


    # Remove NaNs if they exist

    dt = dt.groupby('SMILES').agg(agg_funcs).reset_index()
    dt.columns = ['_'.join(col).strip() for col in dt.columns.values]
    dt.rename(columns={
        'mass_from_formula_<lambda>': 'mz',
        'trj_end_any': 'trj_end_ever',
        'formula_<lambda>': 'formula',
        'SMILES_': 'SMILES'
    }, inplace=True)

    dt['trj_end_ever'].fillna(False, inplace=True)

    # The RMassBank::dbe function was not translated, as its Python equivalent is not readily available
    # You'll need a suitable function to replace it if necessary
    # dt['DBE'] = ...

    dt.dropna(subset=['SMILES'], inplace=True)
    dt.to_csv("Node_table.csv", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process and aggregate CSV files.")
    parser.add_argument("--input_path", type=str, help="Path to the directory containing the CSV files.")
    parser.add_argument("--edge_output", type=str, default="Edge_table.csv",
                        help="Filename for the edge table output. (default: Edge_table.csv)")
    parser.add_argument("--node_output", type=str, default="Node_table.csv",
                        help="Filename for the node table output. (default: Node_table.csv)")

    args = parser.parse_args()

    # Passing arguments to main function
    main(args.input_path)
