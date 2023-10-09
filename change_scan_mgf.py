#python ~/QCXMS/scripts/change_scan_mgf.py merged_predicted_spectra.mgf 
import argparse
from pyteomics import mgf

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mgf_file1', type=str, help='mgf_file1')


    args = parser.parse_args()

    mgf_file1 = args.mgf_file1
  

    add_scan_number = 432

    with open(mgf_file1, 'r') as f1, open('merged_predicted_spectra_changed.mgf', 'w') as output_file:
        # Loop over each line in file2 and write it to the output file with updated scan numbers
        for line in f1:
            if "SCAN=" in line:
                scan_number = int(line.split("SCAN=")[1].split()[0])
                if scan_number >= 1:
                    new_scan_number = scan_number + add_scan_number
                    line = line.replace("SCAN=" + str(scan_number), "SCAN=" + str(new_scan_number))
            output_file.write(line)






if __name__ == "__main__":
    main()
