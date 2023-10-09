#python ~/QCXMS/scripts/random_bin/count_peaks.py merged_predicted_spectra.mgf
import argparse
from pyteomics import mgf
import os


def parse_mgf(data):
    # Split the data into different ion sections
    ion_sections = data.strip().split('BEGIN IONS')

    parsed_data = []

    for section in ion_sections:
        if not section:
            continue

        ion_data = {
            'pepmass': None,
            'charge': None,
            'scan': None,
            'title': None,
            'mslevel': None,
            'peaks': [],
            'peak_count': 0
        }

        lines = section.split('\n')
        for line in lines:
            line = line.strip()
            if line.startswith('PEPMASS='):
                ion_data['pepmass'] = float(line.split('=')[1])
            elif line.startswith('CHARGE='):
                ion_data['charge'] = line.split('=')[1]
            elif line.startswith('SCAN='):
                ion_data['scan'] = line.split('=')[1]
            elif line.startswith('TITLE='):
                ion_data['title'] = line.split('=')[1]
            elif line.startswith('MSLEVEL='):
                ion_data['mslevel'] = int(line.split('=')[1])
            elif line and '=' not in line and "END IONS" not in line:
                # print("check line: ", line)
                mz, intensity = map(float, line.split())
                ion_data['peaks'].append((mz, intensity))
                ion_data['peak_count'] += 1

        parsed_data.append(ion_data)

    return parsed_data




def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mgf_file1', type=str, help='mgf_file1')

    args = parser.parse_args()

    mgf_file1 = args.mgf_file1

    with open(mgf_file1) as mgf_file:
        data = mgf_file.read()

    parsed_data = parse_mgf(data)
    
    # Printing the parsed data along with peak count for each ion section
    nums_peaks = []
    for ion in parsed_data:
        print(f"Title: {ion['title']}")
        print(f"Number of peaks: {ion['peak_count']}")
        print("----------------------------")
        nums_peaks.append(float(ion['peak_count']))
    print(nums_peaks)






if __name__ == "__main__":
    main()
