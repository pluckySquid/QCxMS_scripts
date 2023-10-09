# python ~/QCXMS/scripts/GNPS_json_to_mgf.py ALL_GNPS.json delete_precursor_17
import argparse
import pandas as pd
import json
from rdkit import Chem
import re

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_json', type=str, help='Path to the output json file')
    parser.add_argument('output_file', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_json = args.input_json
    output_file = args.output_file


    with open(input_json, 'r') as f:
        data = json.load(f)

    # Loop through all the strings in the JSON data
    smile_dict = {}
    smile_mgf_dict = {}
    for item in data:
        #print(item)
        if 'Smiles' in str(item):
            smile = item["Smiles"]
            Adduct = item["Adduct"]
            # set the index for the smile
            # if smile in smile_dict.keys():

            if smile == " " or smile == 'N/A':
                continue
            
            if "M+H" not in Adduct:
                print(Adduct)
            else:
                if smile in smile_dict.keys():                   
                    #print("smile ", smile_dict[smile])
                    smile_dict[smile] = smile_dict[smile] + 1
                else:
                    smile_dict[smile] = 1
                
                peak_intensity = item["peaks_json"]
                precursor = item["Precursor_MZ"]

                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, peak_intensity)

                # Convert the matches to float and group them into pairs
                coords = [[float(matches[i]), float(matches[i+1])] for i in range(0, len(matches), 2)]

                mz = [x[0] for x in coords]
                #print(mz)
                intensity = [x[1] for x in coords]

                out_mz_intensity = []
                for i in range(len(mz)):
                    if mz[i] < float(precursor) - 17:
                        #print("float(precursor): ", float(precursor), "mz: ", mz[i])
                        out_mz_intensity.append([mz[i], intensity[i]])


                smile_mgf_dict[str(smile) + "_" + str(smile_dict[smile])] = str(out_mz_intensity) + " id= " + item["SpectrumID"]

    smile_peak_atom_dict = {}
    for key, value in smile_mgf_dict.items():
        #print(key, value)
        if "N/A" not in str(key) and str(key) != "":
            # str(key).split("_")[0] is the smile before adding _#
            try:
                mol = Chem.MolFromSmiles(str(key).split("_")[0])
                mol_withH = Chem.AddHs(mol)
                atom_number = mol_withH.GetNumAtoms()
                smile_peak_atom_dict[str(key)+ ": " + value] = atom_number
            except:
                print("wrong mol", str(key))
            # Replace 'file_path.txt' with the name and path of the file you want to create
    smile_peak_atom_dict = dict(sorted(smile_peak_atom_dict.items(), key=lambda item: item[1], reverse=False))

    #print(smile_peak_atom_dict)

    
    #print dict
    with open(output_file, 'w') as f:
        f.write('smiles: peaks_json, atom\n')

        for key, value in smile_peak_atom_dict.items():
            #print(key, value)
            if "N/A" not in str(key):
                f.write(key + ": " + str(value) + "\n")



if __name__ == "__main__":
    main()