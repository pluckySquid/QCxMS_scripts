import argparse
from fragmentation_py import FragmentEngine
from pyteomics import mgf
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_mz(fragment_mass, charge):
    return fragment_mass / charge

def smiles_to_molblock(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToMolBlock(mol)
    return None

def calculate_precursor_mass(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    # Check if the SMILES can be converted to a molecule
    if mol is None:
        return None
    
    # Calculate the molecular weight (precursor mass)
    mass = Descriptors.ExactMolWt(mol)
    return mass
    
def write_mgf(file, fragment_sets, precursor_list, charge=1):
    print("file", file)
    counter = 1
    with open(file, 'w') as mgf_file:
        for index, fragment_engine in enumerate(fragment_sets):
            mgf_file.write("BEGIN IONS\n")
            
            mgf_file.write(f"PEPMASS={precursor_list[counter]}\n")
            mgf_file.write(f"CHARGE={charge}+\n")
            mgf_file.write(f"SCAN={counter}\n")
            mgf_file.write(f"TITLE=fragment_prediction={counter}\n")
            mgf_file.write(f"MSLEVEL=2\n")
            #mgf_file.write(f"RTINSECONDS=0.0\n")
            #mgf_file.write(f"SEQ={fragment_smiles}\n")  # Include SMILES in SEQ field
            #mgf_file.write("TAXONOMY=unknown\n")
            #mgf_file.write(f"SCORE={score}\n")
            #mgf_file.write(f"COMP={bondbreaks}\n")
            
            mz_set = []
            for fragment_info in fragment_engine.fragment_info:    
                fragment, score, bondbreaks = fragment_info
                title = f"Fragment_{index}_{fragment}"

                fragment_smiles = fragment_engine.get_fragment_info(fragment, 0)[2]
                fragment_mass = fragment_engine.calc_fragment_mass(fragment)

                # Calculate m/z for the fragment
                mz = calculate_mz(fragment_mass, charge)
                mz_set.append(mz)

            mz_set.sort()
            #mgf_file.write(f"totalpeaks = {len(mz_set)}\n")
            for mz in mz_set:
                mgf_file.write(f"{mz} 1.0\n")  # Example intensity value, you can adjust it
            mgf_file.write("END IONS\n\n")

            counter += 1

def predict_spectra(smiles_list, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge):
    fragment_sets = []
    for smiles in smiles_list:
        print("smile", smiles)
        molblock = smiles_to_molblock(smiles)
        # Create a FragmentEngine object for each SMILES string
        fragment_engine = FragmentEngine(molblock, max_broken_bonds, max_water_losses, ionisation_mode, skip_fragmentation, molcharge)

        # Check if the object was accepted
        if not fragment_engine.accepted():
            print(f"Failed to create the FragmentEngine object for SMILES: {smiles}")
            continue

        # Generate fragments
        fragment_count = fragment_engine.generate_fragments()
        print(f"Generated {fragment_count} fragments for SMILES: {smiles}")

        # Get all fragment information
        fragment_sets.append(fragment_engine)

    return fragment_sets

def main(args):
    input_file = args.input_file
    output_mgf = args.output_mgf

    smiles_list = []
    precursor_list = [0]

    # Parse the input file (assuming it's a text file with one SMILES string per line)
    with open(input_file, 'r') as text_file:
        for line in text_file:
            line = line.strip()
            if "smile" not in line:
                smiles_list.append(line.split(",")[0])
                precursor_list.append(calculate_precursor_mass(line.split(",")[0]))
    print("smiles_list", smiles_list)
    print("precursor_list:", precursor_list)

    if not smiles_list:
        print("No SMILES found in the input file.")
        return

    # Predict spectra for each SMILES
    fragment_sets = predict_spectra(smiles_list, args.max_broken_bonds, args.max_water_losses, args.ionisation_mode, args.skip_fragmentation, args.molcharge)

    if not fragment_sets:
        print("No valid spectra predicted.")
        return

    # Write the spectra to an MGF file
    write_mgf(output_mgf, fragment_sets, precursor_list, args.charge)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict spectra for SMILES strings and write to MGF file.")
    parser.add_argument("input_file", type=str, help="Input file containing SMILES strings, one per line")
    parser.add_argument("output_mgf", type=str, help="Output MGF file to write the spectra")

    parser.add_argument("--max_broken_bonds", type=int, default=2, help="Maximum number of broken bonds")
    parser.add_argument("--max_water_losses", type=int, default=1, help="Maximum number of water losses")
    parser.add_argument("--ionisation_mode", type=int, default=1, help="Ionisation mode")
    parser.add_argument("--skip_fragmentation", action="store_true", help="Skip fragmentation")
    parser.add_argument("--molcharge", type=int, default=0, help="Molecule charge")
    parser.add_argument("--charge", type=int, default=1, help="Charge for the spectra")

    args = parser.parse_args()
    main(args)
