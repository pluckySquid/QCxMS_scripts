#python ~/QCXMS/scripts/change_title_mgf.py merged_predicted_spectra.mgf
import argparse
from pyteomics import mgf
import os

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mgf_file1', type=str, help='mgf_file1')

    args = parser.parse_args()

    mgf_file1 = args.mgf_file1

    with open(mgf_file1, 'r') as f1, open('titled_predicted_spectra_changed.mgf', 'w') as output_file:
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

            if "TITLE=QCXMS_SCAN=" in line:
                line = line.replace("TITLE=QCXMS_SCAN=" + str(scan_number), "TITLE=QCXMS_SCAN=" + str(file_contents))
            output_file.write(line)






if __name__ == "__main__":
    main()
