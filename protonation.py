from rdkit import Chem
import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mol_smile', type=str, help='input mole smile')

    args = parser.parse_args()
    
    mol_smile = args.mol_smile

    smiles = []  # Initialize an empty list to store SMILES strings

    with open(mol_smile, 'r') as f:
        lines = f.readlines()

        for line in lines:
            if "smile" not in line:  # Check if "smile" (case insensitive) is not in the line
                smiles.append(line.strip())
    print(smiles)
    
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

    print(protonated_smiles)

if __name__ == "__main__":
    main()