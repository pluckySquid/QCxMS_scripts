#python ~/QCXMS/scripts/random_bin/generate_random_peaks.py merged_predicted_spectra.mgf
import argparse
from pyteomics import mgf
import os
import random



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
                #ion_data['title'] = line.split('=')[1:]
                extracted_str = line.split('=')[1:]
                extracted_str = "=".join(extracted_str)
                ion_data['title'] = extracted_str
            elif line.startswith('MSLEVEL='):
                ion_data['mslevel'] = int(line.split('=')[1])
            elif line and '=' not in line and "END IONS" not in line:
                # print("check line: ", line)
                mz, intensity = map(float, line.split())
                if mz > 50 and mz < ion_data['pepmass'] - 17 and intensity > 0.1:
                    ion_data['peaks'].append((mz, intensity))
                    ion_data['peak_count'] += 1

        parsed_data.append(ion_data)

    return parsed_data


def generate_random_peaks(num_peaks=10, max_peak = 200):
    peaks = []
    for _ in range(num_peaks):
        mz = random.uniform(50, max_peak)  # Random m/z value between 50 and 2000
        intensity = random.uniform(1, 1)  # Random intensity between 0 and 1000
        peaks.append((mz, intensity))
    
    peaks.sort(key=lambda x: x[0])
    return peaks

def create_mgf_with_random_peaks(parsed_data, nums_peaks):
    new_mgf_data = ""
    counter = 0 
    for ion in parsed_data:
        new_mgf_data += "BEGIN IONS\n"
        new_mgf_data += f"PEPMASS={ion['pepmass']}\n"
        new_mgf_data += f"CHARGE={ion['charge']}\n"
        new_mgf_data += f"SCAN={ion['scan']}\n"
        new_mgf_data += f"TITLE={ion['title']}\n"
        new_mgf_data += f"MSLEVEL={ion['mslevel']}\n"

        random_peaks = generate_random_peaks(nums_peaks[counter], ion['pepmass']-17)
        for mz, intensity in random_peaks:
            new_mgf_data += f"{mz:.4f} {intensity:.4f}\n"

        new_mgf_data += "END IONS\n\n"

        counter += 1
    return new_mgf_data

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
        nums_peaks.append(int(ion['peak_count']))
    print(nums_peaks)


    # Write to a new file
    new_mgf_data = create_mgf_with_random_peaks(parsed_data, nums_peaks)
    with open("generated_random_peaks.mgf", "w") as file:
        file.write(new_mgf_data)





if __name__ == "__main__":
    main()
