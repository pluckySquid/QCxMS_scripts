# python ~/QCXMS/scripts/clean_input_dict.py ~/QCXMS/workDir/neccessary_files/FDA-Part2/delete_precursor_20 mols_first_10_mols.txt cleaned_delete_precursor_file
import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_data_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('input_txt_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('output_file', type=str, help='output_json_smile_file is the txt file that has the information about which all the inpt smiles')

    args = parser.parse_args()

    input_data_file = args.input_data_file
    input_txt_file = args.input_txt_file
    output_file = args.output_file

    smile_list = []
    with open(input_txt_file, "r") as input_file:
         for line in input_file:
              if "smiles,energy,trajectories,atom" not in line:
                   smile_list.append(line.split(",")[0])
    print(smile_list)

    with open(input_data_file, "r") as input_file:
        with open(output_file, "w") as out_file:
            for line in input_file:
                    if str(line) != "smiles: peaks_json, atom\n":
                        smile = line.split(":")[0].split("_")[0]
                        if smile in smile_list:
                                print("find smile in smile list:" , smile)
                                out_file.write(line)
                     
                    
if __name__ == "__main__":
    main()
