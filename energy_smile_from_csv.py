#python ~/QCXMS/scripts/energy_smile_from_csv.py ALL_GNPS_cleaned.csv smile_cd.txt
import csv
import argparse
from rdkit import Chem

instrument_list = ["qTof", "Hybrid FT", "Q-TOF", "qToF", "Q-Exactive",
           "Orbitrap", "qTOF", "Maxis II HD Q-TOF Bruker", "Q-Exactive Plus",
           "Q-Exactive Plus Orbitrap Res 14k", "Q-Exactive Plus Orbitrap Res 70k",
           "Maxis HD qTOF", "LC-ESI-QTOF", "LC-ESI-ITFT",
           "LC-ESI-ITTOF", "LC-ESI-QFT", "HPLC-ESI-TOF",
           "ESI-ITFT", "LC-Q-TOF/MS", "ESI-FTICR", "UPLC-ESI-QTOF"]

def get_id_energy_dict(real_energy_file):
    smile_energy_dict = {}
    smile_id_dict = {}

    with open(real_energy_file, "r") as file:
        reader = csv.DictReader(file)
        header = next(reader)
        print(header)
        for row in reader:
            smile = row['Smiles']  # Access the value in the "Smiles" column
            cd = row['collision_energy']
            adduct = row['Adduct']
            GNPS_Inst = row['GNPS_Inst']
            if cd != '' and adduct == 'M+H':
                if not any(GNPS_Inst.lower() == string.lower() for string in instrument_list):
                    continue
                if smile in smile_energy_dict.keys():
                    smile_energy_dict[smile].append(cd)
                    smile_id_dict[smile].append(row['spectrum_id'])
                else:
                    smile_energy_dict[smile] = [cd]
                    smile_id_dict[smile] = [row['spectrum_id']]
    return smile_energy_dict, smile_id_dict


def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_csv_file', type=str, help='input_mgf_file is the mgf file from HPCC')
    parser.add_argument('output_file', type=str, help='input_mgf_file is the mgf file from HPCC')

    args = parser.parse_args()

    input_csv_file = args.input_csv_file
    output_file = args.output_file

    smile_energy_dict, smile_id_dict = get_id_energy_dict(input_csv_file)

    with open(output_file, "w") as out_file:
        for smile in smile_energy_dict:
            mol = Chem.MolFromSmiles(str(smile))
            atom_number = 100
            try: 
                    mol_withH = Chem.AddHs(mol)
                    atom_number = mol_withH.GetNumAtoms()
            except:
                print ("wrong smile")
                continue

            if len(smile_energy_dict[smile]) >= 3 and atom_number <= 30:
                out_file.write(smile + str(smile_energy_dict[smile])+ str(atom_number) +'\n')

    with open("id_" + output_file, "w") as out_file:
        for smile in smile_id_dict:
            mol = Chem.MolFromSmiles(str(smile))
            atom_number = 100
            try: 
                    mol_withH = Chem.AddHs(mol)
                    atom_number = mol_withH.GetNumAtoms()
            except:
                print ("wrong smile")
                continue

            if len(smile_id_dict[smile]) >= 3 and atom_number <= 30:
                out_file.write(smile + str(smile_id_dict[smile])+ str(atom_number) +'\n')
    #print("smile_energy_dict: ", smile_energy_dict)

if __name__ == "__main__":
    main()