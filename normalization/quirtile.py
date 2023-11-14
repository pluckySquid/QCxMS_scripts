import argparse
import numpy as np
import math

# 40 eV:
min_peaks_list_40 = [0, 20, 112, 47, 0, 6, 125, 173, 116, 0, 0, 0, 0, 75, 13, 24, 0, 30, 137, 53]

min_peaks_list_30 = [0, 18, 73, 14, 0, 6, 22, 116, 93, 0, 0, 0, 0, 0, 0, 0, 0, 0, 103, 11]

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
                #print(f"Molecule {molecule['SCAN']}: Number of mz/intensity pairs = {num_peaks}")
            elif line:
                mz, intensity = map(float, line.split())
                if mz > 50 and mz < molecule["PEPMASS"] - 17:
                    molecule["peaks"].append((mz, intensity))

    correct_order_melecules = []

    for i in range(1, len(molecules) + 1):
        for mol in molecules:
            if mol["SCAN"] == i:
                correct_order_melecules.append(mol)
                continue
        
    #print("correct_order_melecules", correct_order_melecules)
    
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
    quantile_list = []


    for i in range(0, len(molecules), group_size):
        group = molecules[i:i + group_size]
        min_peaks = result[molecule_number][1]
        print("min_peaks: ", min_peaks)


        #print("min_peaks: ", min_peaks)
        #print("result[molecule_number]: ", result[molecule_number])
        #print("min_peaks_list_40[molecule_number]: ", min_peaks_list_40[molecule_number])
        quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1]
        # round up
        quantiles_as_int = [math.ceil(min_peaks * q) for q in quantiles]
        print("check quantiles_as_int", quantiles_as_int)
        #quantiles_as_int = [min(math.ceil(min_peaks * q), min_peaks_list_40[molecule_number], min_peaks_list_50[molecule_number], min_peaks_list_60[molecule_number]) for q in quantiles]
        # round down
        #quantiles_as_int = [int(min_peaks * q) for q in quantiles]

        print("quantiles_as_int: ", quantiles_as_int)
        quantile_list.append(quantiles_as_int[1])

        metadata = group[0]  # Take metadata from the first molecule in the group
        result.append((metadata, min_peaks))
        # Output the processed molecules for each molecule in the group
        
        for quantile_index in range(0, 5):
            print("quantile_index: ", quantile_index)
            for molecule in group:
                #print("molecule:", molecule)
                with open(output_file+"_"+str(quantile_index)+".mgf", "a") as processed_file:
                    processed_file.write(f"BEGIN IONS\n")
                    for key, value in molecule.items():
                        if key != "peaks":
                            processed_file.write(f"{key}={value}\n")
                    
                    # Choose the top min_peaks intensity values and corresponding mz values
                    sorted_peaks = sorted(molecule["peaks"], key=lambda x: x[1], reverse=True)
                    if molecule["SCAN"] % group_size >= 3 or molecule["SCAN"] % group_size == 0: # 60 eV +
                        print("molecule[SCAN] % 8 == ?: ", i)
                        top_peaks = sorted_peaks[quantiles_as_int[quantile_index]:quantiles_as_int[quantile_index+1]]
                    else:
                        if molecule["SCAN"] % group_size == 1:
                            print("molecule[SCAN] % 8 == 1: ", i)
                            quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1]
                            quantiles_as_int_temp = [math.ceil(min_peaks_list_30[int(i/group_size)] * q) for q in quantiles]
                            top_peaks = sorted_peaks[quantiles_as_int_temp[quantile_index]:quantiles_as_int_temp[quantile_index+1]]
                        if molecule["SCAN"] % group_size == 2:
                            print("molecule[SCAN] % 8 == 2: ", i)
                            quantiles = [0, 0.2, 0.4, 0.6, 0.8, 1]
                            quantiles_as_int_temp = [math.ceil(min_peaks_list_40[int(i/group_size)] * q) for q in quantiles]
                            top_peaks = sorted_peaks[quantiles_as_int_temp[quantile_index]:quantiles_as_int_temp[quantile_index+1]]
                    # Sort by mz when writing to the file                   
                    top_peaks = sorted(top_peaks, key=lambda x: x[0])
                    #print("quantiles_as_int", quantiles_as_int)
                    #print("quantiles_as_int", quantiles_as_int[quantile_index+1])
                    print("quantiles_as_int[quantile_index]:quantiles_as_int[quantile_index+1]: ", quantiles_as_int[quantile_index], quantiles_as_int[quantile_index+1])
                    print("toppeaks:", len(top_peaks))

                    processed_file.write(f"peaks={len(top_peaks)}\n")
                    
                    for mz, intensity in top_peaks:
                        processed_file.write(f"{mz:.6f} {intensity:.6f}\n")
                    
                    processed_file.write("END IONS\n\n")

        print("min_peaks: ", min_peaks)
        molecule_number += 1
    print("quantile_list: ", quantile_list)
    for i in quantile_list:
        print(i)




def main():
    parser = argparse.ArgumentParser(description="Count mz/intensity pairs in mass spectrometry data")
    parser.add_argument("input_file", help="Input file containing mass spectrometry data")
    parser.add_argument("group_size", type=int, help="Number of molecules per group")
    parser.add_argument("output_file", help="Output file to save results")
    args = parser.parse_args()

    # QCxMS prediction smallest
    # ecom = 6
#     min_peaks_list = [99, 78, 11, 77, 91, 54, 109, 48, 65, 128, 143, 118, 95, 81, 191, 102, 124, 91, 40, 78]
# # ecom 200 traj
#     min_peaks_list = [160, 129, 17, 114, 137, 77, 168, 54, 77, 185, 207, 186, 129, 115, 273, 177, 195, 122, 53, 141]
# # energy 200 traj
#     min_peaks_list = [160, 129, 14, 92, 137, 77, 145, 54, 71, 185, 207, 186, 129, 115, 273, 165, 186, 119, 51, 135]
# # correct:
# # energy test: 50 eV
#     min_peaks_list = [6, 15, 11, 17, 2, 5, 37, 23, 47, 1, 1, 5, 5, 18, 2, 3, 1, 11, 20, 20]
#     # energy_test: 60 eV:
#     min_peaks_list = [18, 40, 12, 30, 3, 14, 76, 24, 72, 13, 3, 14, 11, 34, 5, 7, 7, 51, 21, 44]
#     # 40 eV:
#     #min_peaks_list = [0, 2, 11, 4, 0, 2, 17, 23, 24, 0, 0, 0, 0, 6, 1, 1, 0, 3, 16, 5]
#     # 50eV: 
#     min_peaks_list = [6, 15, 11, 17, 2, 5, 37, 23, 47, 1, 1, 5, 5, 18, 2, 3, 1, 11, 20, 20]

#     # 70 eV:
#     min_peaks_list = [50, 88, 18, 82, 20, 23, 96, 31, 74, 33, 9, 34, 28, 48, 22, 21, 17, 98, 23, 76]

# 50eV, no filter
    min_peaks_list = [88, 143, 115, 114, 21, 35, 224, 182, 258, 14, 20, 40, 72, 164, 18, 34, 29, 128, 146, 203]
    molecules = parse_mz_intensity_pairs(args.input_file)
    # for i in molecules:
    #     print(i)
    results = process_molecules(molecules, args.group_size, min_peaks_list)
    save_results(molecules, args.group_size, args.output_file, results)

if __name__ == "__main__":
    main()
