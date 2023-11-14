import argparse

n_time_crest = [233, 99, 33, 122, 121, 95, 195, 85, 94, 295, 297, 249, 185, 126, 364, 258, 260, 220, 61, 216]

def extract_top_intensities(input_file, output_file):
    # Read the input MGF file
    with open(input_file, "r") as f:
        lines = f.readlines()

    # Initialize a list to store the data lines for each spectrum
    data_lines = []

    # Initialize a variable to keep track of the data lines
    inside_data_block = False

    # Initialize a variable to count the spectra
    spectrum_count = 0

    # Initialize variables to store the "BEGIN IONS" and "END IONS" lines
    begin_ions_line = None
    end_ions_line = None

    # Initialize a variable to store the "TITLE" line
    title_line = None

    # Initialize a variable to store the "PEPMASS" line
    pepmass_line = None

    # Initialize variables to store the "CHARGE," "SCAN," and "MSLEVEL" lines
    additional_info = []

    # Iterate through the lines and extract the data
    molecule_num = 0
    for line in lines:
        line = line.strip()
        if line == "BEGIN IONS":
            inside_data_block = True
            spectrum_count += 1
            begin_ions_line = line
        elif line == "END IONS":
            inside_data_block = False
            end_ions_line = line
        elif inside_data_block:
            # Store the "TITLE" line
            if line.startswith("TITLE="):
                title_line = line
            # Store the "PEPMASS" line
            elif line.startswith("PEPMASS="):
                pepmass_line = line
            # Store the "CHARGE," "SCAN," and "MSLEVEL" lines
            elif line.startswith(("CHARGE=", "SCAN=", "MSLEVEL=")):
                additional_info.append(line)
            else:
                data_lines.append(line)

        # If we have processed the data for one spectrum, sort and write the top 200 intensities to the output file
        if not inside_data_block and data_lines:
            # Sort the data lines by intensity, only if the lines have the expected format
            valid_data_lines = [line for line in data_lines if len(line.split()) == 2]
            valid_data_lines.sort(key=lambda x: float(x.split()[1]), reverse=True)
            
            # Keep only the top 200 intensities
            valid_data_lines = valid_data_lines[:n_time_crest[molecule_num]]
            valid_data_lines = sorted(valid_data_lines, key=lambda x: float(x.split()[0]))


            # Write the "BEGIN IONS" line
            with open(output_file, "a") as f:
                f.write(begin_ions_line + "\n")

                # Write the "PEPMASS" line
                f.write(pepmass_line + "\n")

                # Write the "TITLE" line
                f.write(title_line + "\n")

                # Write the "CHARGE," "SCAN," and "MSLEVEL" lines
                for line in additional_info:
                    f.write(line + "\n")

                # Write the selected data lines
                for line in valid_data_lines:
                    f.write(line + "\n")

                # Write the "END IONS" line
                f.write(end_ions_line + "\n\n")

            # Clear the data_lines list and additional_info for the next spectrum
            data_lines = []
            additional_info = []
            print(n_time_crest[molecule_num])
            molecule_num += 1

    print(f"Top 200 intensities for {spectrum_count} spectra written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract top intensities from an MGF file.")
    parser.add_argument("input_file", help="Input MGF file")
    parser.add_argument("output_file", help="Output MGF file")

    args = parser.parse_args()
    
    # Clear the output file before writing to it
    with open(args.output_file, "w") as f:
        f.write("")  # Clear the file

    extract_top_intensities(args.input_file, args.output_file)
