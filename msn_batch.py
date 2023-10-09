import random
import sys

def main(input_file):
    # Read all lines from the input file
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Remove the first line (header) from the list
    header = lines.pop(0)

    # Shuffle the lines randomly
    random.shuffle(lines)

    # Split the shuffled lines into three groups
    group1 = lines[:100]
    group2 = lines[100:300]
    group3 = lines[300:]

    # Write the groups to their respective output files
    output_file_1 = 'batch1.txt'
    output_file_2 = 'batch2.txt'
    output_file_3 = 'batch3.txt'

    with open(output_file_1, 'w') as outfile:
        outfile.write(header)
        outfile.writelines(group1)

    with open(output_file_2, 'w') as outfile:
        outfile.write(header)
        outfile.writelines(group2)

    with open(output_file_3, 'w') as outfile:
        outfile.write(header)
        outfile.writelines(group3)

    print("Files generated successfully.")

if __name__ == "__main__":
    # Check for the correct number of command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py input.txt")
        sys.exit(1)

    # Get the input file name from the command-line argument
    input_file = sys.argv[1]

    # Call the main function with the input file
    main(input_file)
