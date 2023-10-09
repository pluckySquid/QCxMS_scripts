import csv

def check_if_gnps_library(line):
    items = line.split(',')
    return items[-2] == "GNPS-SELLECKCHEM-FDA-PART2"

input_file = 'ALL_GNPS_cleaned.csv'
output_file = 'small.csv'

matching_lines = []

with open(input_file, 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if check_if_gnps_library(','.join(row)):
            matching_lines.append(row)

with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(matching_lines)
