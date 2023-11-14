import os
import tarfile
import subprocess
import shutil

def xyz_to_smile(xyz):
    subprocess.check_call(['obabel', xyz,  '-O','%s' %(xyz + ".smi")])
    with open(xyz+".smi", 'r') as file:
        content = file.read()
    return content


def untar_file(tar_file, destination_folder):
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(destination_folder)
    print("done un tar process")

def convert_to_smiles(input_folder, output_file):
    subfolder_path = os.path.join(input_folder, 'TMPQCXMS')
    
    folders = os.listdir(subfolder_path)
    with open(output_file, 'w') as output:
        for folder in folders:
            folder_path = os.path.join(subfolder_path, folder)
            files = os.listdir(folder_path)
            
            for file in files:
                if file.endswith('.xyz') and file.count('.') == 3:
                    # Check if the file matches the corrected specified pattern
                    print("file:", file)
                    input_file = os.path.join(folder_path, file)
                    print("file:", input_file)
                    smi = xyz_to_smile(input_file)
                    print(input_file)
                    smile = smi.split()[0]
                    print("smile:", smile)
                    output.write(f"{file}\t{smile}\n")
                    # Remove the generated .smi files
                    smi_files = [f for f in os.listdir(folder_path) if f.endswith('.smi')]
                    for smi_file in smi_files:
                        os.remove(os.path.join(folder_path, smi_file))
def main():
    # Find the current working directory and join with 'ecom test'
    main_folder = os.path.join(os.getcwd(), 'ecom_test')

    # Path to the folder containing tar files
    tar_folder = os.path.join(main_folder, 'pqcxms_grouped')

    # Path to the folder where the SMILES files will be stored
    output_folder = 'output_smiles'
    os.makedirs(output_folder, exist_ok=True)

    for tar_file in os.listdir(tar_folder):
        if tar_file.endswith('.tar'):
            untar_folder = os.path.join(main_folder, f'TMPQCXMS_{tar_file}')  # Create a unique folder for each tar file

            untar_file(os.path.join(tar_folder, tar_file), untar_folder)
            
            # Extract trajectory number from the tar file name
            traj_number = tar_file.split('_')[-1].split('.')[0]

            print("untar_folder: ", untar_folder)

            convert_to_smiles(untar_folder, os.path.join(output_folder, f'{traj_number}_smiles.txt'))

            # Clean up: Remove the untarred folder
            shutil.rmtree(untar_folder)

    print("Conversion completed successfully.")

if __name__ == "__main__":
    main()
