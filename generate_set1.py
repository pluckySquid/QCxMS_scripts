#python ~/QCXMS/scripts/generate_set1.py mols_set1.txt mols_multi-or-crest_tbd.txt 
import argparse
import pandas as pd
import json
from rdkit import Chem
def cal_proteomer_numers(mol):
    # Find heteroatoms (N,O,S) in molecule
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in ['N', 'O', 'S']]

    # Protonate each heteroatom and generate new SMILES
    protonated_smiles = []
    for atom in heteroatoms:
        atom.SetFormalCharge(1)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        protonated_smiles.append(Chem.MolToSmiles(mol, isomericSmiles=True))
        # reset the atom to the original state
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
    return len(protonated_smiles)

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_txt', type=str, help='Path to the output json file')
    parser.add_argument('output_txt', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_txt = args.input_txt
    output_txt = args.output_txt

    smiles_list=[]
    with open(input_txt, 'r') as f:
        data = f.read().splitlines()
    
        for line in data:
            if "smiles" not in line:
                smiles_list.append(line.split(",")[0])

    #print("smiles_list:", smiles_list)

    total_time = 0
    mol_atom_dict = {}
    for key in smiles_list:
        
        if "N/A" not in str(key) and str(key) != ' ':
            mol = Chem.MolFromSmiles(str(key))
            try: 
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                #if atom_number >= 36 and atom_number <= 40:
                if atom_number <= 50:
                    mol_atom_dict[str(key)] = atom_number
                    mol_time =  8 * (0.16 * atom_number * atom_number - 2.77 * atom_number)
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
    with open("mols_n_times_crest.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5," + str(100 * cal_proteomer_numers(mol)) + ",n2,ecom," + str(atom_number) + "\n"  

                f.write(string)

    total_time = 0
    with open("mols_multi_proteomers.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                total_time += cal_proteomer_numers(mol) * (0.16 * atom_number * atom_number - 2.77 * atom_number)
                f.write(string)
        print("mols_multi_proteomers total_time: ", total_time)

    with open("mols_ecom_test.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",6,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",6.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",7,100,n2,ecom," + str(atom_number) + "\n"

                f.write(string)


    with open("mols_energy_test.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",30,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",40,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",50,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",60,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",70,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",80,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",90,100,n2,elab," + str(atom_number) + "\n"
                string = string + str(key) + ",100,100,n2,elab," + str(atom_number) + "\n"

                f.write(string)
                # Replace 'file_path.txt' with the name and path of the file you want to create
        #f.write("total_atom_number" + str(total_atom_number))
    with open("mols_reproducibility_20_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,20,n2,ecom," + str(atom_number) + "\n"
                f.write(string)

    with open("mols_reproducibility_40_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,40,n2,ecom," + str(atom_number) + "\n"
                f.write(string)

    with open("mols_reproducibility_5.5_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,5.5,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

    with open("mols_reproducibility_80_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,80,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

    with open("mols_reproducibility_100_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,100,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

    with open("mols_reproducibility_200_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,200,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

    with open("mols_reproducibility_300_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,300,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

    with open("mols_reproducibility_500_traj.txt", 'w') as f:
        f.write('smiles,energy,trajectories,gas,energy_type,atom\n')

        for key, value in sorted_mol_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                mol = Chem.MolFromSmiles(str(key))
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                string = str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"
                string = string + str(key) + ",5.5,500,n2,ecom," + str(atom_number) + "\n"

                f.write(string)

if __name__ == "__main__":
    main()
