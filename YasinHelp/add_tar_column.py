import csv
import argparse
import multiprocessing
from rdkit import Chem
from rdkit.Chem import AllChem

def protonation(smiles):
    # Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

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

def add_tar_info_row(row, molecule_index, process_type):
    tar_file = f"{molecule_index}-{1}_pqcxms.tar"  # Assuming the second-to-last column is the gas type
    row.append(tar_file)
    return row

def add_tar_info_row_multi(row, molecule_index, process_type):
    smiles = row[0]
    proteomer_num = protonation(smiles)
    tar_file = ""
    for i in range(1, proteomer_num+1):
        tar_file += f"{molecule_index}-{i}_pqcxms.tar; "  
        
    row.append(tar_file)
    return row

def add_tar_info(input_file, output_file, process_type):
    with open(input_file, 'r') as input_csv, open(output_file, 'w', newline='') as output_csv:
        reader = csv.reader(input_csv, delimiter=',')
        writer = csv.writer(output_csv, delimiter=',')
        
        # Add header with the new column
        header = next(reader)
        header.append('tar_file')
        writer.writerow(header)
        
        # Optional multiprocessing
        if process_type == 'multi':
            pool = multiprocessing.Pool()
            molecule_index = 0
            for row in pool.starmap(add_tar_info_row_multi, [(row, molecule_index, process_type) for molecule_index, row in enumerate(reader)]):
                writer.writerow(row)
                molecule_index += 1
            pool.close()
            pool.join()
        else:
            molecule_index = 0
            for row in reader:
                row = add_tar_info_row(row, molecule_index, process_type)
                writer.writerow(row)
                molecule_index += 1

def main():
    parser = argparse.ArgumentParser(description='Add tar file information to a CSV file.')
    parser.add_argument('input_file', help='Input CSV file')
    parser.add_argument('output_file', help='Output CSV file with tar file information')
    parser.add_argument('--process_type', default='single', choices=['single', 'multi'], help='Processing type (single or multi, default is single)')

    args = parser.parse_args()

    add_tar_info(args.input_file, args.output_file, args.process_type)

if __name__ == "__main__":
    main()
