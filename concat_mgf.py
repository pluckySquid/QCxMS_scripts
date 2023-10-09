import argparse
from pyteomics import mgf

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mgf_file1', type=str, help='mgf_file1')
    parser.add_argument('mgf_file2', type=str, help='mgf_file2')


    args = parser.parse_args()

    mgf_file1 = args.mgf_file1
    mgf_file2 = args.mgf_file2 

    max_scan_number = 1
    with open(mgf_file1, 'r') as f1:
        for line in f1:
            if "SCAN=" in line:
                # Extract the scan number from the line
                scan_number = int(line.split("SCAN=")[1].split()[0])
                # Update the maximum scan number if necessary
                if scan_number > max_scan_number:
                    max_scan_number = scan_number
    # Read the contents of the first file
    max_scan_number = max_scan_number
    print(max_scan_number)

    with open(mgf_file1, 'r') as f1, open(mgf_file2, 'r') as f2, open('merged_predicted_spectra.mgf', 'w') as output_file:
        # Loop over each line in file1 and write it to the output file
        for line in f1:
            output_file.write(line)

        # Loop over each line in file2 and write it to the output file with updated scan numbers
        for line in f2:
            if "SCAN=" in line:
                scan_number = int(line.split("SCAN=")[1].split()[0])
                new_scan_number = max_scan_number + scan_number # + 8 # if the last molecule at the end is a wrong molecule
                line = line.replace("SCAN=" + str(scan_number), "SCAN=" + str(new_scan_number))
            output_file.write(line)


if __name__ == "__main__":
    main()
