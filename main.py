import sys
import argparse
from rdkit import Chem, rdBase
from xyz2mol import *
import re
import os
import numpy as np

def parse_out_file(filename):
    parsed_data = {}
    current_key = None
    reading_charges = False
    reading_fragments = False
    current_line = None

    with open(filename, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            last_line = current_line
            current_line = lines[i].strip()
            previous_line = lines[i - 1].strip() if i - 1 > 0 else None
            next_line = lines[i+1].strip() if i+1 < len(lines) else None

            if current_line:  # only process line if it's not empty
                if "Heating  trajectory" in current_line:
                    current_key = "heating_trajectory" + str(
                        len([k for k in parsed_data.keys() if k.startswith('heating_trajectory')]) + 1)
                    parsed_data[current_key] = {"charges_per_fragment": [], "fragments": {}, "step": []}
                elif "trajectory" in current_line and "collision" in current_line:
                    current_key = "CID" + str(len([k for k in parsed_data.keys() if k.startswith('CID')]) + 1)
                    parsed_data[current_key] = {"charges_per_fragment": [], "fragments": {}, "step": []}
                elif "- Entering Mean-Free-Path simulation -" in current_line:
                    current_key = "MD" + str(len([k for k in parsed_data.keys() if k.startswith('MD')]) + 1)
                    parsed_data[current_key] = {"charges_per_fragment": [], "fragments": {}, "step": []}

                if "FRAGMENTATION occured!" in current_line:
                    a = 0
                    while not previous_line or not previous_line.split()[0].isdigit():
                        a -= 1
                        previous_line = lines[i+a].strip() if i+a > 0 else None
                    step_start = int(previous_line.split()[0]) if previous_line else None


                    a = 0
                    while not next_line or not next_line.split()[0].isdigit():
                        a += 1
                        next_line = lines[i+a].strip() if i+a < len(lines) else None
                    step_end = int(next_line.split()[0]) if next_line else None


                    parsed_data[current_key]['step'].append([step_start, step_end])

                if "Summed charges per fragment" in current_line:
                    reading_charges = True
                    reading_fragments = False
                elif "mass                formula" in current_line:
                    reading_charges = False
                    reading_fragments = True
                elif reading_charges:
                    data = current_line.split()
                    fragment_id, charge = data[0], data[-1]
                    # Check if the first part of the line can be converted to an integer.
                    # If not, we've reached the end of the "Summed charges per fragment" section.
                    try:
                        fragment_id = int(fragment_id)
                        charge = float(charge)
                    except ValueError:
                        reading_charges = False
                        continue

                    parsed_data[current_key]["charges_per_fragment"].append({"id": fragment_id, "charge": charge})
                elif reading_fragments:
                    data = current_line.split()
                    if current_line.startswith('M='):
                        keys = ["mass", "formula", "q", "pop", "spin", "|q IPB|", "diss time (ps)"]
                        data[0] = data[0][2:]
                        for i, key in enumerate(keys):
                            if key not in parsed_data[current_key]["fragments"]:
                                parsed_data[current_key]["fragments"][key] = [data[i]]
                            else:
                                parsed_data[current_key]["fragments"][key].append(data[i])

    return parsed_data

def add_matrix(matrices, matrix, num_matrices):
    matrices.append(matrix)
    if len(matrices) > num_matrices:
        matrices.pop(0)
    return matrices


def separateXYZGroupIntoXYZIndividualStrings(xyzGroupedMolecules):
    firstLine = xyzGroupedMolecules[0]
    firstLineNoSpaces = firstLine.strip()
    num_atoms = int(firstLineNoSpaces)
    numLinesPerXYZ = num_atoms + 2
    numXYZStructures = int(len(xyzGroupedMolecules) / numLinesPerXYZ)
    xyzStrings = []
    for startXYZ in range (0,numXYZStructures):
        xyzString=""
        for indexLineFile in range((startXYZ*numLinesPerXYZ),((startXYZ+1)*numLinesPerXYZ)):
            xyzString = xyzString + xyzGroupedMolecules[indexLineFile] + '\n'
        xyzStrings.append(xyzString)
    return xyzStrings


def writeInAFile(fileName, content):
    f = open(fileName, "w")
    f.writelines(content)
    f.close()

def appendInAFile(fileName, content):
    f = open(fileName, "a")
    f.writelines(content+'\n')
    f.close()

def longest_substring(input_string):
    # Split the string into a list by the '.' character
    substrings = input_string.split('.')

    # Filter the list for substrings that contain '+'
    substrings = [s for s in substrings if '+' in s]

    # If there are no substrings with '+', return an empty string
    if not substrings:
        return ''
    # Find the longest substring that contains '+'
    longest_substring = max(substrings, key=len)

    return longest_substring


def get_Smiles_with_charge(smiles, charge = 1):
    # Split the string into individual molecules
    molecules = smiles.split('.')
    # Calculate and print the total formal charge for each molecule
    for mol_smiles in molecules:
        rdBase.DisableLog('rdApp.warning')
        mol = Chem.MolFromSmiles(mol_smiles)
        rdBase.EnableLog('rdApp.warning')
        total_charge = Chem.rdmolops.GetFormalCharge(mol)
        if total_charge == charge:
            return mol_smiles
    return('')

def main(xyz_path, csv_path):


    xyz_filename = xyz_path.split('/')[-1]
    CID_substring = re.compile(r'CID\d+')
    is_CID_file = CID_substring.search(xyz_filename) != None
    num_matrices = 4

    with open(xyz_path) as fp:

        contentXYZFile = fp.read().splitlines()
        listOfXYZStrings = separateXYZGroupIntoXYZIndividualStrings(contentXYZFile)
        i =0
        prev_smiles = None
        fileName = "temp.xyz"
        fileName = os.path.join(os.path.dirname(csv_path), fileName)
        fileNameForStructures = csv_path
        AC_prev = []
        dMat_prev = []

        for xyz_string in listOfXYZStrings:
            failed = ''
            writeInAFile(fileName,xyz_string)

            try:
                atoms, charge, xyz_coordinates = read_xyz_file(fileName)
            except Exception as e:
                failed = 'xyz'
                print(str(xyz_path) + ' was not read successfully!')
                print(e)


            try:
                if failed == '':
                    molStructure, AC, dMat = xyz2mol(atoms, xyz_coordinates, charge=1, allow_charged_fragments=True,
                                use_graph=True, use_huckel=False, embed_chiral=False,
                                use_atom_maps=False, tr_previous_AC = AC_prev, tr_previous_dMat = dMat_prev, N2collision = is_CID_file)

                    AC_prev = add_matrix(AC_prev, AC, num_matrices)
                    dMat_prev = add_matrix(dMat_prev, dMat, num_matrices)
                    if molStructure != 'same':
                        smiles = Chem.MolToSmiles(molStructure[0])
                        smiles = harmonize_smiles_rdkit(smiles)


            except Exception as e:
                failed = 'mol'
                print(str(i) + 'molecule processing failed! (' + xyz_path + ')')
                print(e)

            if failed != '':
                appendInAFile(fileNameForStructures, str(i) + ',' + failed)
                print(i, failed)
            elif i == 0 or smiles != prev_smiles:

                appendInAFile(fileNameForStructures, str(i) + ',' + smiles)
                prev_smiles = smiles
                print(i, smiles)

            i += 1


if __name__ == "__main__":
    #import cProfile
    #profiler = cProfile.Profile()
    #profiler.enable()

    parser = argparse.ArgumentParser(description='Process XYZ and write to CSV')
    parser.add_argument('--xyz_path', '-xyz', type=str, required=True, help='Path to XYZ file')
    parser.add_argument('--csv_path', '-csv', type=str, required=True, help='Path to CSV file for output')
    args = parser.parse_args()

    main(args.xyz_path, args.csv_path)

    #main("C:/PostDoc/Ming_time/example_files/water.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #main("C:/PostDoc/Ming_time/example_files/single_proton.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #main("C:/PostDoc/Ming_time/example_files/CID3.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #main("C:/PostDoc/Ming_time/example_files/CID4.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #print('next')
    #main("C:/PostDoc/Ming_time/example_files/MDtrj.27.4.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #main("C:/PostDoc/Ming_time/example_files/MDtrj.75.2.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")
    #main("C:/PostDoc/Ming_time/example_files/MDtrj.30.3.xyz", "C:/PostDoc/Ming_time/example_files/53_CID1_test.csv")

    #profiler.disable()
    #profiler.print_stats(sort="time")
