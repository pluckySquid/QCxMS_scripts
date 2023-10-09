#python ~/QCXMS/scripts/sort_json.py ../../GNPS-SELLECKCHEM-FDA-PART2.json mols_32_atoms.txt
import argparse
import pandas as pd
import json
from rdkit import Chem

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_json', type=str, help='Path to the output json file')
    parser.add_argument('output_txt', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_json = args.input_json
    output_txt = args.output_txt

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
    for item in data:
        #print(item)
        if 'smiles' in str(item):
            #print("str(item): ", str(item))
            smile = item["smiles"]
            print("smile: ", smile)
            #print("item[molecularFormula]: ",item["molecularFormula"])
            #print("item[name]: ",item["name"])
            continue
        if 'headline' in str(item):
            Adduct = item["headline"]
            print(Adduct)
            # if smile in smile_dict.keys():
            if smile in smile_dict.keys() and "[M+H]+" in Adduct:    
                print("smile ", smile_dict[smile])
                smile_dict[smile] = smile_dict[smile] + 1
            elif "M+H" in Adduct:
                smile_dict[smile] = 1
            else:
                continue

            if smile_dict[smile] > 1:
                print(smile)
    

    total_time = 0
    mol_atom_dict = {}
    for key, value in smile_dict.items():
        #print(key, value)
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



    #print dict
    total_atom_number = 0
    with open(output_txt, 'w') as f:
        f.write('smiles,energy,trajectories,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",50,100," + str(atom_number) + "\n"
                string = string + str(key) + ",60,100," + str(atom_number) + "\n"
                string = string + str(key) + ",70,100," + str(atom_number)  + "\n"
                string = string + str(key) + ",80,100," + str(atom_number)  + "\n"
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
