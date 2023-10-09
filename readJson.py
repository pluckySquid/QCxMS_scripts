import argparse
import pandas as pd
import json


def main():
    parser = argparse.ArgumentParser(description='convert the csv result to json file')
    parser.add_argument('input_json', type=str, help='Path to the output json file')
    parser.add_argument('output_txt', type=str, help='index for the molecule')

    args = parser.parse_args()

    input_json = args.input_json
    output_txt = args.output_txt


    with open(input_json, 'r') as f:
        data = json.load(f)

    # Loop through all the strings in the JSON data
    smile_dict = {}
    for item in data:
        #print(item)
        if 'Smiles' in str(item):
            smile = item["Smiles"]
            if smile in smile_dict.keys():
                #print("smile ", smile_dict[smile])
                smile_dict[smile] = smile_dict[smile] + 1
            else:
                smile_dict[smile] = 1
    
    #print dict
    with open(output_txt, 'w') as f:
        f.write('smiles,energy\n')

        for key, value in smile_dict.items():
            print(key, value)
            if "N/A" not in str(key):
                string = str(key) + ",10\n"
                string = string + str(key) + ",20\n"
                string = string + str(key) + ",30\n"
                string = string + str(key) + ",40\n"
                string = string + str(key) + ",50\n"
                string = string + str(key) + ",60\n"
                string = string + str(key) + ",70\n"
                f.write(string)
                # Replace 'file_path.txt' with the name and path of the file you want to create
        




if __name__ == "__main__":
    main()