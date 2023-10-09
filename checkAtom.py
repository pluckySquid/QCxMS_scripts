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

with open('small_mol.txt', 'r') as file:
    lines = file.readlines()
    for row in lines:

        #print(row[1])
        smile = row
        #print(smile)
        #smile = "CCCCCCn1cccc1C(=O)c1cccc2ccccc12"
        
        # Create an RDKit molecule object from the SMILES string
        mol = Chem.MolFromSmiles(smile)
        mol_withH = Chem.AddHs(mol)

        # Get the atom of interest, e.g., the first carbon atom
        #atom = mol.GetAtomWithIdx(0)

        # Get the atomic number of the atom
        atom_number = mol_withH.GetNumAtoms()
        print(atom_number)







