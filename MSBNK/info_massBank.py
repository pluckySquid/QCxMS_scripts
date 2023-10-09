#python ~/QCXMS/scripts/sort_json.py ../../GNPS-SELLECKCHEM-FDA-PART2.json mols_32_atoms.txt
import argparse
import pandas as pd
import json
from rdkit import Chem
import collections
import os

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["data_count", "headline", "monoisotopicMolecularWeight", "identifier", "spectrum"]
)

instrument_list = ["ESI-ITTOF", "ESI-ITFT", "ESI-QTOF", "ESI-TOF", "LC-ESI-ITFT", "LC-ESI-ITTOF", "LC-ESI-QFT", "LC-ESI-QTOF", "LC-ESI-TOF"]

folder_path = '/home/user/QCXMS/data/massbank/spectra/MassBank-MassBank-data-67b1750'

def find_spec(identifier, monoisotopicMolecularWeight, delete_precursor_switch):
    folder_name = identifier.split("-")[1]
    print("folder_name: ", folder_name)
    file_name = identifier + str(".txt")
    file_path = os.path.join(folder_path, folder_name,file_name)
    mz = []
    intensity = []
    real_intensity = []
    mz_intensity = []

    if os.path.isfile(file_path):
        #print("File exists!")
        annotation_flag = False
        with open (file_path) as file:
            for line in file:
                if "PK$NUM_PEAK:" in line:
                    annotation_flag = True
                    num_peaks= line.split("PK$NUM_PEAK: ")[1]
                    #print("found PK$NUM_PEAK:")
                    continue
                if "//\n" == line:
                    #print("end of file")
                    break
                if annotation_flag and "PK$PEAK:" not in line:
                    #print("line.split()", line.split())
                    if not delete_precursor_switch:
                        mz.append(line.split()[0])
                        intensity.append(line.split()[1])
                        real_intensity.append(line.split()[2])
                        mz_intensity.append((line.split()[0], line.split()[2]))
                    if delete_precursor_switch and float(line.split()[0]) < monoisotopicMolecularWeight + 1 - 20:
                        mz.append(line.split()[0])
                        intensity.append(line.split()[1])
                        real_intensity.append(line.split()[2])
                        mz_intensity.append((line.split()[0], line.split()[2]))

            #print("mz ", mz)
            #print("intensity ", intensity)
            #print("real_intensity ", real_intensity)
            # if not annotation_flag:
            #     print("no annotation ", file_path)

    else:
        print("File does not exist.", file_path)
    return mz_intensity

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_json', type=str, help='Path to the output json file')
    parser.add_argument('output_txt', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_json = args.input_json
    output_txt = args.output_txt

    delete_precursor_switch = True
    count = 1
    with open(input_json, 'r') as f:
        for line in f:
            #print(line)
            count = count  + 1 
            if count > 1000:
                break

    with open(input_json, 'r') as f:
        data = json.load(f)

    # Loop through all the strings in the JSON data
    smile_dict = {}
    counter = 0
    smile = ""
    monoisotopicMolecularWeight = ""
    identifier = ""
    Headline = ""
    for item in data:
        #print(item)
        if 'smiles' in str(item):
            #print("str(item): ", str(item))
            smile = item["smiles"]
            monoisotopicMolecularWeight = item["monoisotopicMolecularWeight"]
            identifier = item["identifier"]
            spectrum = find_spec(identifier,monoisotopicMolecularWeight,delete_precursor_switch)
            #print(Headline)
            # if smile in smile_dict.keys():
            if smile in smile_dict.keys() and "[M+H]+" in Headline:    
                #print("smile ", smile_dict[smile])
                smile_dict[smile].append(SpectrumTuple(len(smile_dict[smile]) + 1, Headline, monoisotopicMolecularWeight, identifier, spectrum))
            elif "[M+H]+" in Headline:
                smile_dict[smile] = [SpectrumTuple(1,Headline, monoisotopicMolecularWeight, identifier, spectrum)]
            else:
                continue

            if len(smile_dict[smile]) > 1:
                print("more than 1 smile: ", smile)
            print("smile: ", smile)
            print("spectrum: ", spectrum)
            #print("item[molecularFormula]: ",item["molecularFormula"])
            #print("item[name]: ",item["name"])
            continue
        if 'headline' in str(item):
            if any(instrument in item["headline"] for instrument in instrument_list):
                #print("instrument: ", item["headline"])
                headline = "||. "
                Headline = headline + item["headline"] + " .||"
                

                
    

    total_time = 0
    mol_atom_dict = {}
    for key, value in smile_dict.items():
        #print("key, value", key, value)
        if "N/A" not in str(key) and str(key) != ' ':
            mol = Chem.MolFromSmiles(str(key))
            try: 
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                #if atom_number >= 36 and atom_number <= 40:
                if atom_number <= 40:
                    mol_atom_dict[str(key)] = atom_number
                    mol_time =  0.04269 * atom_number ** 2 - 0.09106 * atom_number - 5.28
                    total_time = total_time + mol_time
            except: 
                with open("wrong_mols.txt", "w") as wrong_mols_file:
                    wrong_mols_file.write(str(key))
            # Replace 'file_path.txt' with the name and path of the file you want to create
    sorted_mol_atom_dict = dict(sorted(mol_atom_dict.items(), key=lambda item: item[1], reverse=False))
    print("total_time: ", total_time)
    #print(sorted_mol_atom_dict)

    sorted_smile_headline_dict = {}
    for smile in sorted_mol_atom_dict.keys():
        sorted_smile_headline_dict[smile] = smile_dict[smile]


    #print dict
    total_atom_number = 0
    with open(output_txt, 'w') as f:
        f.write('smiles,atom,headline\n')

        for key, value in sorted_smile_headline_dict.items():
            print("key, value", key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) +","+ str(atom_number) + "," + str(value) +"\n"

                # string = string + str(key) + ",70,100," + str(atom_number)  + "\n"
                # string = string + str(key) + ",80,100," + str(atom_number)  + "\n"
                # string = string + str(key) + ",90,100," + str(atom_number)  + "\n"
                # string = string + str(key) + ",100,100," + str(atom_number)  + "\n"
                # string = str(key)  + ",10," + str(atom_number) + "\n"
                # string = string + str(key) + ",13," + str(atom_number) + "\n"
                # string = string + str(key) + ",16," + str(atom_number)  + "\n"
                # string = string + str(key) + ",19," + str(atom_number)  + "\n"
                # string = string + str(key) + ",22," + str(atom_number)  + "\n"
                # string = string + str(key) + ",25," + str(atom_number)  + "\n"
                # string = string + str(key) + ",28," + str(atom_number)  + "\n"
                # string = string + str(key) + ",31," + str(atom_number) + "\n"
                # string = string + str(key) + ",34," + str(atom_number)  + "\n"
                # string = string + str(key) + ",37," + str(atom_number)  + "\n"
                # string = string + str(key) + ",40," + str(atom_number)  + "\n"
                total_atom_number += atom_number
                f.write(string)
                # Replace 'file_path.txt' with the name and path of the file you want to create
        f.write("total_atom_number" + str(total_atom_number))

if __name__ == "__main__":
    main()
