import subprocess
import glob
import argparse
import os
import re
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs

def xyz_to_smile(xyz):
    subprocess.check_call(['obabel', xyz,  '-O','%s' %(xyz + ".smi")])
    with open(xyz+".smi", 'r') as file:
        content = file.read()
    return content

def get_mass(smiles):
    try:   
        mol = Chem.MolFromSmiles(smiles)
        mass = Descriptors.ExactMolWt(mol)
        # H_mass = 1.007277
        return mass # + H_mass
    except: 
        return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read smiles and write mol files')
    parser.add_argument('key', type=str, help='Path to the output folder')

    args = parser.parse_args()

    key = args.key

    pattern = r'\d+\.\d+\.\d+\.xyz'
    folder = 'TMPQCXMS'

    # Use glob to search for matching files in all subdirectories
    matching_files = []
    for root, dirs, files in os.walk(folder):
        for filename in files:
            if re.match(pattern, filename):
                matching_files.append(os.path.join(root, filename))

    mass_list = []
    for file in matching_files:
        smi = xyz_to_smile(file)
        print(smi)
        smile = smi.split()[0]
        print("smile:", smile)
        smile_mass = get_mass(smile)
        print("mass:", smile_mass)
        if smile_mass not in mass_list:
            mass_list.append(smile_mass)

    print(mass_list)
    with open("mass_" + key + ".txt", "w") as file:
        file.write(str(mass_list))



