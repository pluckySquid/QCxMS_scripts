import argparse
import pandas as pd
import json
from rdkit import Chem
import re

def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_file', type=str, help='Path to the output json file')
    parser.add_argument('output_file', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file


    with open(input_file, 'r') as file:
        lines = file.readlines()

    
    modified_lines = []

    # Process each line
    line_counter = 1
    for line in lines:
        # Define the regular expression pattern
        pattern = r"SpectrumTuple\(data_count=\d+, headline='(.+?)', monoisotopicMolecularWeight=(\d+\.\d+), identifier='(.+?)', spectrum=\[(.*?)\]\)"

        # Find all occurrences of the pattern in the data
        matches = re.findall(pattern, line)

        print("matches:", matches)

        # Extract the relevant information for each SpectrumTuple
        spectrum_tuples = []
        for match in matches:
            data_count, headline, weight, identifier, spectrum_data = match
            spectrum_values = re.findall(r"\('([\d.]+)', '(\d+)'\)", spectrum_data)
            spectrum_tuples.append({
                "data_count": data_count,
                "headline": headline,
                "monoisotopicMolecularWeight": float(weight),
                "identifier": identifier,
                "spectrum": [(float(mz), int(intensity)) for mz, intensity in spectrum_values if float(mz) < (float(weight) + 1 - 20)] 
            })

        # Print the extracted SpectrumTuples
        for spectrum_tuple in spectrum_tuples:
            print(spectrum_tuple)

        modified_lines.append(spectrum_tuples)
        line_counter += 1
        if line_counter > 4:
            break

    # Save the modified lines to a new file
    with open(output_file, 'w') as file:
        file.writelines(str(modified_lines))


if __name__ == "__main__":
    main()