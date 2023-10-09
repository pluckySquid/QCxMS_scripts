import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('output_file', type=str, help='json_smile_peak_file is the converted file that has the information that each smiles in the json file with its peak and intensity')

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    with open(input_file, "r") as file:
        lines = file.readlines()

    filtered_lines = []
    for line in lines:
        parts = line.split()
        if "BEGIN IONS" in line or "END IONS" in line or "PEPMASS" in line or "CHARGE" in line or "SCAN" in line or "TITLE" in line or "MSLEVEL" in line:
            filtered_lines.append(line)
            continue
        print("parts: ", parts)

        if len(parts) >= 2 and float(parts[1]) > 0.1:
            filtered_lines.append(line)

    with open(output_file, "w") as file:
        file.writelines(filtered_lines)

if __name__ == "__main__":
    main()
