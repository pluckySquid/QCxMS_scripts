#python ~/QCXMS/scripts/csv_to_txt.py csv table.txt
import argparse
import pandas as pd
import json
from rdkit import Chem
import csv

# Sort the SMILES by atom numbers
def count_atoms(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        mol_withH = Chem.AddHs(mol)
        atom_number = mol_withH.GetNumAtoms()
    except:
        print("wrong smile:", smiles)
    if mol is None:
        return 0
    return atom_number

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_csv', type=str, help='Path to the output json file')
    #parser.add_argument('output_txt', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_csv = args.input_csv
   # output_txt = args.output_txt

    #Read the CSV file and extract SMILES
    smiles_list = []
    with open(input_csv, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if len(row) > 0 and  "SMILES" not in row:
                smiles = row[0]
                smiles_list.append(smiles)

    sorted_smiles_list = sorted(smiles_list, key=count_atoms)

    # Write each sorted SMILES to a separate text file
    
    with open(f'sorted_MSn.txt', 'w') as txt_file:
        txt_file.write("smiles,energy,trajectories,atom")
        for i, smiles in enumerate(sorted_smiles_list):
            #print(key, value)
            if "N/A" not in str(smiles):
                mol = Chem.MolFromSmiles(str(smiles))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(smiles) + ",80,100," + str(atom_number)+ "\n"
                #total_atom_number += atom_number
                txt_file.write(string)





if __name__ == "__main__":
    main()
