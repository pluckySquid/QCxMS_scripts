# python ~/QCXMS/scripts/all_in_one.py ~/QCXMS/workDir/neccessary_files/GNPS-SELLECKCHEM-FDA-PART2.json ~/QCXMS/workDir/neccessary_files/real_energy.txt ../GNPS-SELLECKCHEM-FDA-PART2_less_than_40_atoms 2
import argparse
import pandas as pd
import json
from rdkit import Chem
import collections

SmileTuple = collections.namedtuple(
    "smile", ["SpectrumID", "atomNumber", "realEnergy", "real_PeakIntensity"]
)

def write_file(sorted_dict, n):
    min_atom_numbers = {}
    with open("output.txt", "w") as f:    
        for key, value in sorted_dict:
            if value.realEnergy not in min_atom_numbers:
                min_atom_numbers[value.realEnergy] = [(value.atomNumber, key, value)]
            else:
                # current_min_atom_number, _ , _  = min_atom_numbers[value.realEnergy]
                # if value.atomNumber < current_min_atom_number:
                #     print("wrong sort")
                min_atom_numbers[value.realEnergy].append((value.atomNumber, key, value))

        for energy, list_of_tuple in min_atom_numbers.items():
            print("len(list_of_tuple), energy", len(list_of_tuple), energy)
            if len(list_of_tuple) >= n:
                (min_atom_number, key, value) = list_of_tuple[n-1]
            else: 
                print("not enough", energy)
                (min_atom_number, key, value) = list_of_tuple[0]
            f.write(f"SpectrumID: {value.SpectrumID}, atomNumber: {min_atom_number}, realEnergy: {energy}, real_PeakIntensity: {value.real_PeakIntensity}, smile: {key}\n")

    total_time = 0
    with open("mols_experiment.txt", "w") as f:
        f.write('smiles,energy,atom,realEnergy\n')
        for energy, list_of_tuple in min_atom_numbers.items():  
            if len(list_of_tuple) >= n:
                (min_atom_number, key, value) = list_of_tuple[n-1]
            else: 
                print("not enough experiments", energy)
                (min_atom_number, key, value) = list_of_tuple[0]
            if min_atom_number <= 70: 
                string = str(key)  + ",10," + str(min_atom_number) + "," + energy + "\n"
                string = string + str(key) + ",20," + str(min_atom_number) + "," + energy + "\n"
                string = string + str(key) + ",30," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",40," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",50," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",60," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",70," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",80," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",90," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",100," + str(min_atom_number)  + "," + energy + "\n"
                string = string + str(key) + ",110," + str(min_atom_number)  + "," + energy + "\n"     
                f.write(string)
                mol_time =  0.04269 * min_atom_number ** 2 - 0.09106 * min_atom_number - 5.28
                total_time = total_time + mol_time * 7
    print("total_time", total_time)

# my_dict[scan_number] = [smile, energy, atom_numbers]
def get_energy(real_energy_file, SpectrumID):
    with open(real_energy_file, "r") as file:
        for line in file:
            if line != "	spectrum_id	collision_energy\n":
                id = line.split()[1]
                energy = line.split()[2]
                if id == SpectrumID:
                    return energy


def read_json(input_json, input_energy, output_txt, batch_number):

    with open(input_json, 'r') as f:
        data = json.load(f)

    # Loop through all the strings in the JSON data
    smile_dict = {}
    smile_info_dict = {}
    total_time = 0
    for item in data:
        if 'Smiles' in str(item):
            smile = item["Smiles"]
            Adduct = item["Adduct"]
            if smile in smile_dict.keys() and "M+H" in Adduct:    
                print("smile ", smile_dict[smile])
                smile_dict[smile] = smile_dict[smile] + 1
            elif "M+H" in Adduct:
                smile_dict[smile] = 1
                if "N/A" not in smile and smile != ' ':
                    mol = Chem.MolFromSmiles(smile)
                    try: 
                        mol_withH = Chem.AddHs(mol)
                        atom_number = mol_withH.GetNumAtoms()
                        smile_info_dict[smile] = SmileTuple(item["SpectrumID"], atom_number, get_energy(input_energy, item["SpectrumID"]), item["peaks_json"])
                        mol_time =  0.04269 * atom_number ** 2 - 0.09106 * atom_number - 5.28
                        if atom_number <= 40:
                            total_time = total_time + mol_time
                    except: 
                        with open("wrong_mols.txt", "w") as wrong_mols_file:
                            wrong_mols_file.write(smile)        
            else:
                continue

            if smile_dict[smile] > 1:
                print("multiple ", smile)

    
    sorted_smile_info_dict = dict(sorted(smile_info_dict.items(), key=lambda item: (item[1].realEnergy, item[1].atomNumber), reverse=False))
    print("total_time: ", total_time)
    #print(sorted_smile_info_dict)

    write_file(sorted_smile_info_dict.items(), batch_number)

    # total_atom_number = 0
    # with open("mols_experiment.txt", 'w') as f:
    #     f.write('smiles,energy,atom\n')
    #     for key, value in sorted_smile_info_dict.items():
    #         if "N/A" not in str(key):
    #             mol = Chem.MolFromSmiles(str(key))
    #             # string = str(key)  + ",10," + str(atom_number) + "\n"
    #             # string = string + str(key) + ",20," + str(atom_number) + "\n"
    #             # string = string + str(key) + ",30," + str(atom_number)  + "\n"
    #             # string = string + str(key) + ",40," + str(atom_number)  + "\n"
    #             # string = string + str(key) + ",50," + str(atom_number)  + "\n"
    #             # string = string + str(key) + ",60," + str(atom_number)  + "\n"
    #             # string = string + str(key) + ",70," + str(atom_number)  + "\n"
    #             string = str(key)  + ",10," + str(value.atomNumber) + "energe: "+ value.realEnergy + "\n"
    #             total_atom_number += atom_number
    #             f.write(string)

    return sorted_smile_info_dict
        
def get_scan_smile_dict(input_json_smile_file): 
    scan_smile_dict = {}
    with open(input_json_smile_file) as input_json_smile_file:
        scanNumber = 0
        for line in input_json_smile_file:
            if line != "smiles,energy,atom\n":
                scan_smile_dict[str(scanNumber)] = line.split(",")[0]
            scanNumber = scanNumber + 1

    return scan_smile_dict

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_json_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('input_energy_file', type=str, help='json_smile_peak_file is the converted file that has the information that each smiles in the json file with its peak and intensity')
    parser.add_argument('input_experiment_file', type=str, help='input_mgf_file is the mgf file from HPCC')
    parser.add_argument('batch_number', type=str, help='input_mgf_file is the mgf file from HPCC')

    args = parser.parse_args()

    input_json_file = args.input_json_file
    input_energy_file = args.input_energy_file 
    input_experiment_file = args.input_experiment_file
    batch_number = int(args.batch_number)

    # input_mgf_file = args.input_mgf_file
    # real_energy_file = args.real_energy_file

    sorted_smile_info_dict = read_json(input_json_file, input_energy_file, "mols_experiment.txt", batch_number)

    scan_smile_dict = get_scan_smile_dict(input_experiment_file)
    #print(scan_smile_dict)




if __name__ == "__main__":
    main()
