#python ~/QCXMS/scripts/get_table.py merged_predicted_spectra.mgf
import argparse
from pyteomics import mgf
import os
import csv


def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mols_txt', type=str, help='input_mgf')
    parser.add_argument('input_mgf', type=str, help='input_mgf')

    args = parser.parse_args()

    mols_txt = args.mols_txt
    input_mgf = args.input_mgf

    qcxms_in = ""

    #change mgf scan to proteomers first: 
    with open(input_mgf, 'r') as f1, open('titled_predicted_spectra_changed.mgf', 'w') as output_file:
        # Loop over each line in file2 and write it to the output file with updated scan numbers
        for line in f1:
            scan_number = 'mol-protonation'
            file_contents = 'mol_name'
            if "SCAN=" in line:
                scan_number = line.split("SCAN=")[1].split()[0]
                print("scan_num:", scan_number)
                mol_index = int(scan_number.split("-")[0])
                protonation_index = int(scan_number.split("-")[1])

                file_dir = os.getcwd() + "/protonation/"
                file_name = str(mol_index - 1) + "_m_" + str(protonation_index) + ".txt"
                file_path = os.path.join(file_dir, file_name)
                print(file_path)

                # Check if the file exists
                if os.path.exists(file_path):
                    # Open the file
                    with open(file_path, 'r') as file:
                        # Read the file contents
                        file_contents = file.read()
                        
                    # Print the file contents
                    print(file_contents)
                else:
                    print("File not found.")

                file_dir = os.getcwd() + "/protonation/"
                file_name = "qcxms.in"
                qcxms_in_path = os.path.join(file_dir, file_name)

                if os.path.exists(file_path):
                    # Open the file
                    with open(qcxms_in_path, 'r') as file:
                        # Read the file contents
                        qcxms_in = file.read().replace('\n', ',')
                else:
                    print("File not found.")
                

            if "TITLE=QCXMS_SCAN=" in line:
                line = line.replace("TITLE=QCXMS_SCAN=" + str(scan_number), "TITLE=QCXMS_SCAN=" + str(file_contents))
            output_file.write(line)




    smile_id = 0
    with open(mols_txt, 'r') as f1, open("titled_predicted_spectra_changed.mgf", 'r') as mgf_file, open('table.csv', 'w', newline='') as output_file:
        csv_writer = csv.writer(output_file)

        # Write the header row
        header = ["count", "smile_id", "smile", "proteomer_id", "proteomer", "CE", "link", "qcxms.in"]
        csv_writer.writerow(header)
        
        # Loop over each line in file2 and write it to the output file with updated scan numbers
        count = 1
        for line in f1:
            if "smiles,energy," not in line:
                smile = line.split(",")[0]
                CE = line.split(",")[1]
                proteomer_id = 1
                file_dir = os.getcwd() + "/protonation/"
                file_name = str(smile_id) + "_m_" + str(proteomer_id) + ".txt"
                file_path = os.path.join(file_dir, file_name)
                #if the proteomer exist
                while os.path.exists(file_path):
                    print("file_path exist", file_path)
                    
                    with open(file_path, 'r') as file:
                    # Read the file contents
                        file_contents = file.read()
                    proteomer_name = file_contents

                    link = "https://fileserver.wanglab.science/QCXMS/" + os.getcwd().split("QCXMS/")[1] + "/pqcxms_grouped/" + str(smile_id) + "-" + str(proteomer_id) + "_pqcxmsms.tar"

                    data = [str(smile_id) + "-" + str(proteomer_id), smile_id, smile, proteomer_id, proteomer_name, CE, link, qcxms_in]
                    csv_writer.writerow(data)

                    proteomer_id += 1
                    file_name = str(smile_id) + "_m_" + str(proteomer_id) + ".txt"
                    file_path = os.path.join(file_dir, file_name)

                    count += 1

                smile_id += 1
                
            # scan_number = 'mol-protonation'
            # file_contents = 'mol_name'
            # if "SCAN=" in line:
            #     scan_number = line.split("SCAN=")[1].split()[0]
            #     print("scan_num:", scan_number)
            #     mol_index = int(scan_number.split("-")[0])
            #     protonation_index = int(scan_number.split("-")[1])

            #     file_dir = os.getcwd() + "/protonation/"
            #     file_name = str(mol_index - 1) + "_m_" + str(protonation_index) + ".txt"
            #     file_path = os.path.join(file_dir, file_name)
            #     print(file_path)

            #     # Check if the file exists
            #     if os.path.exists(file_path):
            #         # Open the file
            #         with open(file_path, 'r') as file:
            #             # Read the file contents
            #             file_contents = file.read()
                        
            #         # Print the file contents
            #         print(file_contents)
            #     else:
            #         print("File not found.")

            # if "TITLE=QCXMS_SCAN=" in line:
            #     line = line.replace("TITLE=QCXMS_SCAN=" + str(scan_number), "TITLE=QCXMS_SCAN=" + str(file_contents))
            # output_file.write(line)






if __name__ == "__main__":
    main()
