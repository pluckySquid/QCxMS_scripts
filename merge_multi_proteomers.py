import argparse
from collections import defaultdict

def merge_proteomers(input_file, output_file):
    # Define a dictionary to store data for each mz value and molecule
    molecule_data = defaultdict(lambda: defaultdict(float))
    molecule_info = defaultdict(dict)  # Store molecule-specific information
    current_molecule = "1"
    previous_molecule = "0"

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            parts = line.split()
            if line.startswith("SCAN="):
                current_molecule = line.split("=")[1].split("-")[0]
            elif len(parts) == 2:
                #print(line)
                if line.startswith("BEGIN IONS"):
                # this is a new proteomer, it can be a new molecule, can be a new proteomer for the previous molecule
                    pass 
                elif line.startswith("END IONS"):
                    if molecule_info[current_molecule].get('PROTEOMERS_NUMS', 0):
                        molecule_info[current_molecule]['PROTEOMERS_NUMS'] += 1
                    else:
                        molecule_info[current_molecule]['PROTEOMERS_NUMS'] = 1
                    print("molecule_info[current_molecule]", current_molecule, molecule_info[current_molecule])
                else:
                    mz, intensity = map(float, parts)
                    if current_molecule:
                        molecule_data[current_molecule][mz] += intensity
            elif line.startswith("PEPMASS"):
                molecule_info[current_molecule]["PEPMASS"] = line.split("=", 1)[1]
            elif line.startswith("CHARGE"):
                molecule_info[current_molecule]["CHARGE"] = line.split("=", 1)[1]
            elif line.startswith("TITLE"):
                molecule_info[current_molecule]["TITLE"] = line.split("=", 1)[1]
            elif line.startswith("MSLEVEL"):
                molecule_info[current_molecule]["MSLEVEL"] = line.split("=", 1)[1]
            # elif line.startswith("END "):
            #     print(line)
            #     if molecule_info[current_molecule]["PROTEOMERS_NUMS"]:
            #         molecule_info[current_molecule]["PROTEOMERS_NUMS"] += 1
            #     else:
            #         molecule_info[current_molecule]["PROTEOMERS_NUMS"] = 1

    # Average the intensities for merged data
    with open(output_file, 'w') as output:
        count = 1
        for molecule, mz_values in molecule_data.items():
            print(count)
            output.write("BEGIN IONS\n")
            output.write(f"PEPMASS={molecule_info[molecule]['PEPMASS']}\n")
            output.write(f"CHARGE={molecule_info[molecule]['CHARGE']}\n")
            output.write(f"SCAN={molecule}\n")
            output.write(f"TITLE={molecule_info[molecule]['TITLE']}\n")
            output.write(f"MSLEVEL={molecule_info[molecule]['MSLEVEL']}\n")
            for mz, intensity in sorted(mz_values.items()):
                avg_intensity = intensity / molecule_info[current_molecule]["PROTEOMERS_NUMS"]
                output.write(f"{mz} {avg_intensity}\n")
            output.write("END IONS\n")
            output.write("\n")  # Add an extra line after END IONS
            count += 1

    print("len(molecule_info)", len(molecule_info))
    for i in molecule_info:
        print("molecule_info[i]", i, molecule_info[i])
        print("molecule_info[i]", molecule_info[i]['PROTEOMERS_NUMS'])

def main():
    parser = argparse.ArgumentParser(description='Merge proteomers and calculate average intensities.')
    parser.add_argument('input_file', help='Input file containing mz, intensity, and proteomer data.')
    parser.add_argument('output_file', help='Output file to save the merged average intensities.')
    args = parser.parse_args()

    merge_proteomers(args.input_file, args.output_file)

if __name__ == '__main__':
    main()
