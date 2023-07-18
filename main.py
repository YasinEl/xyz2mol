import sys

from rdkit import Chem, rdBase
from xyz2mol import *

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
    with open(xyz_path) as fp:
        contentXYZFile = fp.read().splitlines()
        listOfXYZStrings = separateXYZGroupIntoXYZIndividualStrings(contentXYZFile)
        i =0
        prev_smiles = None
        fileName = "temp.xyz"
        fileNameForStructures = csv_path
        AC = []

        for xyz_string in listOfXYZStrings:

            writeInAFile(fileName,xyz_string)

            atoms, charge, xyz_coordinates = read_xyz_file(fileName)

            molStructure, AC = xyz2mol(atoms, xyz_coordinates, charge=1, allow_charged_fragments=True,
                        use_graph=True, use_huckel=False, embed_chiral=False,
                        use_atom_maps=False, tr_previous_AC = AC)

            smiles = Chem.MolToSmiles(molStructure[0])

            #charged_smiles = get_Smiles_with_charge(smiles)



            # Apply the function to a test SMILES string
            #corrected_smiles = correct_charges(smiles)
            #print('cor: ' + corrected_smiles)

            if i == 0 or smiles != prev_smiles:
                appendInAFile(fileNameForStructures, str(i) + ',' + smiles)

            prev_smiles = smiles
            i+=1
            print(i,smiles)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process XYZ and write to CSV')
    parser.add_argument('--xyz_path', type=str, required=True, help='Path to XYZ file')
    parser.add_argument('--csv_path', type=str, required=True, help='Path to CSV file for output')
    args = parser.parse_args()
    main(args.xyz_path, args.csv_path)
