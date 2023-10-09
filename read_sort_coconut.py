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

    smile_dict = {}
    with open(input_json, 'r') as f:
        for line in f:
            if '"smiles":"' in str(line):
                #print(str(line))
                smile = str(line).split('"smiles":"')[1].split('"')[0]
                #print(smile)
                # if smile in smile_dict.keys():
                if smile in smile_dict.keys():    
                    #print("smile ", smile_dict[smile])
                    smile_dict[smile] = smile_dict[smile] + 1
                else:
                    smile_dict[smile] = 1
    

    total_time = 0
    mol_atom_dict = {}
    for key, value in smile_dict.items():
        #print(key, value)
        if "N/A" not in str(key) and str(key) != ' ':
            mol = Chem.MolFromSmiles(str(key))
            try: 
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                mol_atom_dict[str(key)] = atom_number
                mol_time =  0.04269 * atom_number ** 2 - 0.09106 * atom_number - 5.28
                if atom_number <= 70:
                    total_time = total_time + mol_time
            except: 
                with open("wrong_mols_" + str(output_txt) + ".txt", "a") as wrong_mols_file:
                    wrong_mols_file.write(str(key) + "\n")
            # Replace 'file_path.txt' with the name and path of the file you want to create
    sorted_mol_atom_dict = dict(sorted(mol_atom_dict.items(), key=lambda item: item[1], reverse=False))
    print("total_time: ", total_time)
    #print(sorted_mol_atom_dict)



    #print dict
    total_atom_number = 0
    with open(output_txt, 'w') as f:
        f.write('smiles,energy,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key)  + ",10," + str(atom_number) + "\n"
                # string = string + str(key) + ",20," + str(atom_number) + "\n"
                # string = string + str(key) + ",30," + str(atom_number)  + "\n"
                # string = string + str(key) + ",40," + str(atom_number)  + "\n"
                # string = string + str(key) + ",50," + str(atom_number)  + "\n"
                # string = string + str(key) + ",60," + str(atom_number)  + "\n"
                # string = string + str(key) + ",70," + str(atom_number)  + "\n"
                total_atom_number += atom_number
                f.write(string)
                # Replace 'file_path.txt' with the name and path of the file you want to create
        f.write("total_atom_number" + str(total_atom_number))
        




if __name__ == "__main__":
    main()
