import rdkit
import csv
from rdkit import Chem
from rdkit.Chem import Draw
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np



counter = 0
small_counter = 0
small_mol = []

with open('summary.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if counter == 0:
            counter = counter + 1
            continue
        #print(row[1])
        smile = row[3]
        #smile = "CCCCCCn1cccc1C(=O)c1cccc2ccccc12"
        
        # Create an RDKit molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smile)
        mol_withH = Chem.AddHs(mol)

        # Get the atom of interest, e.g., the first carbon atom
        #atom = mol.GetAtomWithIdx(0)

        # Get the atomic number of the atom
        atom_number = mol_withH.GetNumAtoms()
        print(atom_number)

        counter = counter + 1


        small_mol.append(atom_number)

        # if atom_number <= 50:
        #     small_counter += 1
        #     small_mol.append(smile)



with open("small_mol.txt", "w") as txt_file:
    for line in small_mol:
        txt_file.write(str(line) + "\n")
print("overall small molecues are: ", small_counter)

