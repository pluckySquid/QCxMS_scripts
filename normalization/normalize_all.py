import argparse

def parse_mz_intensity_pairs(input_file):
    molecules = []
    with open(input_file, 'r') as file:
        molecule = {}
        for line in file:
            line = line.strip()
            if line.startswith("BEGIN IONS"):
                molecule = {"PEPMASS": None, "CHARGE": None, "SCAN": None, "TITLE": None, "peaks": []}
            elif line.startswith("PEPMASS"):
                molecule["PEPMASS"] = float(line.split('=')[1])
            elif line.startswith("CHARGE"):
                molecule["CHARGE"] = line.split('=')[1]
            elif line.startswith("SCAN"):
                molecule["SCAN"] = int(line.split('=')[1])
            elif line.startswith("TITLE"):
                molecule["TITLE"] = line.split('TITLE=')[1]
            elif line.startswith("MSLEVEL"):
                molecule["MSLEVEL"] = int(line.split('=')[1])
            elif line.startswith("END IONS"):
                molecules.append(molecule)
                num_peaks = len(molecule["peaks"])
                print(f"Molecule {molecule['SCAN']}: Number of mz/intensity pairs = {num_peaks}")
            elif line:
                mz, intensity = map(float, line.split())
                #if intensity > 0.1 and mz > 50 and mz < molecule["PEPMASS"] - 17:
                if mz > 50 and mz < molecule["PEPMASS"] - 17:
                    molecule["peaks"].append((mz, intensity))

    correct_order_melecules = []

    for i in range(1, len(molecules) + 1):
        for mol in molecules:
            if mol["SCAN"] == i:
                correct_order_melecules.append(mol)
                continue
    return correct_order_melecules

def process_molecules(molecules, group_size, min_peaks_list):
    results = []
    molecule_number = 1
    for i in range(0, len(molecules), group_size):
        group = molecules[i:i + group_size]
        min_peaks = min_peaks_list[molecule_number - 1]  # Use the value from the list
        metadata = group[0]  # Take metadata from the first molecule in the group
        results.append((molecule_number, min_peaks))
        # Output the processed molecules for each molecule in the group
        for molecule in group:
            print("molecule: ", molecule)
            with open("processed_molecules.txt", "a") as processed_file:
                processed_file.write(f"BEGIN IONS\n")
                for key, value in molecule.items():
                    if key != "peaks":
                        processed_file.write(f"{key}={value}\n")
                
                # Choose the top min_peaks intensity values and corresponding mz values
                sorted_peaks = sorted(molecule["peaks"], key=lambda x: x[1], reverse=True)
                top_peaks = sorted_peaks[:min_peaks]
                
                # Sort by mz when writing to the file
                top_peaks = sorted(top_peaks, key=lambda x: x[0])
                
                for mz, intensity in top_peaks:
                    processed_file.write(f"{mz:.6f} {intensity:.6f}\n")
                
                processed_file.write("END IONS\n\n")
        #print("min_peaks: ", min_peaks)
        molecule_number += 1
    return results

def save_results(molecules, group_size, output_file, result):
    molecule_number = 0
    for i in range(0, len(molecules), group_size):
        group = molecules[i:i + group_size]
        min_peaks = result[molecule_number][1]
        metadata = group[0]  # Take metadata from the first molecule in the group
        result.append((metadata, min_peaks))
        # Output the processed molecules for each molecule in the group
        for molecule in group:
            with open(output_file, "a") as processed_file:
                processed_file.write(f"BEGIN IONS\n")
                for key, value in molecule.items():
                    if key != "peaks":
                        processed_file.write(f"{key}={value}\n")
                
                # Choose the top min_peaks intensity values and corresponding mz values
                sorted_peaks = sorted(molecule["peaks"], key=lambda x: x[1], reverse=True)
                top_peaks = sorted_peaks[:min_peaks]
                
                
                # Sort by mz when writing to the file
                top_peaks = sorted(top_peaks, key=lambda x: x[0])

                processed_file.write(f"peaks={len(top_peaks)}\n")
                
                for mz, intensity in top_peaks:
                    processed_file.write(f"{mz:.6f} {intensity:.6f}\n")
                
                processed_file.write("END IONS\n\n")
        print("min_peaks: ", min_peaks)
        molecule_number += 1




def main():
    parser = argparse.ArgumentParser(description="Count mz/intensity pairs in mass spectrometry data")
    parser.add_argument("input_file", help="Input file containing mass spectrometry data")
    parser.add_argument("group_size", type=int, help="Number of molecules per group")
    parser.add_argument("output_file", help="Output file to save results")
    args = parser.parse_args()

    #100 traj
    # QCxMS prediction smallest
    min_peaks_list = [99, 78, 11, 77, 91, 54, 109, 48, 65, 128, 143, 118, 95, 81, 191, 102, 124, 91, 40, 78]

    min_peaks_list = [18, 29, 23, 23, 5, 7, 45, 37, 52, 3, 4, 8, 15, 33, 4, 7, 6, 26, 30, 41]  


    molecules = parse_mz_intensity_pairs(args.input_file)
    for i in molecules:
        print(i)
    results = process_molecules(molecules, args.group_size, min_peaks_list)
    save_results(molecules, args.group_size, args.output_file, results)

if __name__ == "__main__":
    main()
