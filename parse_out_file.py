import sys
import argparse
from rdkit import Chem, rdBase
from xyz2mol import *
import re
import os
import json
import numpy as np


def get_decimal_places(num):
    """Returns the number of decimal places in a number."""
    s = str(num)
    if '.' in s:
        return len(s) - s.index('.') - 1
    return 0

def round_to_smaller_decimal_places(num1, num2):
    """Rounds both numbers to the smallest number of decimal places between them."""
    num1 = float(num1)
    num2 = float(num2)
    places = min(get_decimal_places(num1), get_decimal_places(num2))
    return round(num1, places), round(num2, places)


def parse_out_file(filename):
    parsed_data = {'status': ''}
    current_key = None
    reading_charges = False
    reading_fragments = False
    current_line = None
    used_keys = []
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                last_line = current_line
                current_line = lines[i].strip()
                previous_line = lines[i - 1].strip() if i - 1 > 0 else None
                next_line = lines[i+1].strip() if i+1 < len(lines) else None

                if current_key not in used_keys and current_key is not None:
                    used_keys.append(current_key)


                if current_line:  # only process line if it's not empty
                    if "Heating  trajectory" in current_line:
                        current_key = "heating_trajectory" + str(
                            len([k for k in parsed_data.keys() if k.startswith('heating_trajectory')]) + 1)
                        parsed_data[current_key] = {"charges_per_fragment": [], "fragments": {}, "step": []}
                    elif "trajectory" in current_line and "collision" in current_line:
                        current_id = str(len([k for k in parsed_data.keys() if k.startswith('CID')]) + 1)
                        current_key = "CID" + current_id
                        parsed_data[current_key] = {"charges_per_fragment": [], "fragments": {}, "step": []}
                    elif "MFP traj." in current_line:
                        current_key = "MD" + current_id + '_' + current_line[-1]
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
                            charge = charge
                        except ValueError:
                            reading_charges = False
                            continue

                        parsed_data[current_key]["charges_per_fragment"].append({"id": fragment_id, "charge": charge, "used": False})
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

                    if 'statistical charge' in current_line and len(used_keys) > 1:
                        previous_key = used_keys[-2]
                        extracted_value = current_line.split()[-1]
                        for item in parsed_data[previous_key]["charges_per_fragment"]:
                            rounded_charge, rounded_extracted = round_to_smaller_decimal_places(item['charge'],
                                                                                                extracted_value)
                            if rounded_charge == rounded_extracted:
                                item['used'] = True


        #assign last nominal charge
        max_index = max(enumerate(parsed_data[current_key]["charges_per_fragment"]), key=lambda x: float(x[1]["charge"]))[0]
        len_index = len(parsed_data[current_key]["charges_per_fragment"])
        if len_index > 1:
            parsed_data[current_key]["charges_per_fragment"][max_index]['used'] = True
        else:
            parsed_data[current_key]["charges_per_fragment"][0]['used'] = True

    except:
        parsed_data['status'] = 'failed'

    return parsed_data


if __name__ == "__main__":
    #import cProfile
    #profiler = cProfile.Profile()
    #profiler.enable()
    #outfile = parse_out_file('C:/PostDoc/Ming_time/example_files/qcxms.out')
    #from pprint import pprint
    #pprint(outfile)


    parser = argparse.ArgumentParser(description='Parse qcxms.out')
    parser.add_argument('--filepath', '-i', type=str, required=True, help='Path to qcxms.out')
    args = parser.parse_args()

    file_path = args.filepath


    if file_path.endswith('qcxms.out'):
        out_dict = parse_out_file(file_path)
        filename = os.path.basename(file_path).split(".out")[0]

        with open(filename + '.json', 'w') as file:
            json.dump(out_dict, file, indent = 4)


