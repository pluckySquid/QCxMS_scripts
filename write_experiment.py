#python ~/QCXMS/scripts/write_experiment.py 10 3 21
import argparse

smile_list_small =["C1=CC=C2C(=C1)C(C(=CC2=O)C)=O", "C1=C(C2=C(C(=C1Cl)O)N=CC=C2)Cl", "C12=C(C=CC(=C1)OC(F)(F)F)N=C(S2)N", "C12OC(C3C4C1=C(C(OC(C(=C(C=3)O)O)=4)=O)C=C(C=2O)O)=O"]
smile_list_big = ["C1=NC(=CC(=C1)OC2=CC(=C(C=C2)NC(NC3=CC=C(C(=C3)C(F)(F)F)Cl)=O)F)C(=O)NC", "C1(=CN=C2C(=C1)C(=CN2)C(=O)C3=C(C(=CC=C3F)NS(=O)(=O)CCC)F)C4=CC=C(C=C4)Cl", "C1(=CC(=C(N=C1)N)O[C@@H](C2=C(C=CC(=C2Cl)F)Cl)C)C3=CN(N=C3)C4CCNCC4", "N1=C(N=C(C=C1NC2=NC=C(S2)C(NC3=C(C=CC=C3C)Cl)=O)N4CCN(CC4)CCO)C"]
smile_list_small = ["C1=CC=C2C(=C1)C(C(=CC2=O)C)=O", "C1=C(C2=C(C(=C1Cl)O)N=CC=C2)Cl", "C12=C(C=CC(=C1)OC(F)(F)F)N=C(S2)N", "N1=CN=C2C(=C1SC3=C(N=CN3C)[N+](=O)[O-])NC=N2", "C12OC(C3C4C1=C(C(OC(C(=C(C=3)O)O)=4)=O)C=C(C=2O)O)=O", "C(C1C(=CC=C(C=1)Cl)O)(=O)NC2C(=CC(=CC=2)[N+](=O)[O-])Cl", "C1(=C(C=C2C(=C1)S(N(C(N2)CCl)C)(=O)=O)Cl)S(=O)(=O)N", "C1=C(C(=CC=C1Cl)C(OCC2C=CSC=2Cl)CN3C=CN=C3)Cl", "C1=C(C=CC(=C1)C(C2=C(Cl)C=C(C=C2Cl)N3N=CC(NC3=O)=O)C#N)Cl"]
#smile_list_small = ["C1=CC=C2C(=C1)C(C(=CC2=O)C)=O", "C12OC(C3C4C1=C(C(OC(C(=C(C=3)O)O)=4)=O)C=C(C=2O)O)=O", "C1(=C(C=CC(=C1)N2C(N(C(C2(C)C)=O)C3=CC=C(C(=C3)C(F)(F)F)C#N)=S)C(NC)=O)F"]
smile_list_small = ["C1=C(C2=C(C(=C1Cl)O)N=CC=C2)Cl", "C1=C(C(=CC=C1Cl)C(OCC2C=CSC=2Cl)CN3C=CN=C3)Cl", "C1=C(C=CC(=C1)C(C2=C(Cl)C=C(C=C2Cl)N3N=CC(NC3=O)=O)C#N)Cl", "C(C(CS(=O)(=O)C1C=CC(=CC=1)F)(C)O)(=O)NC2C=C(C(=CC=2)C#N)C(F)(F)F", "C1(=C(C=CC(=C1)N2C(N(C(C2(C)C)=O)C3=CC=C(C(=C3)C(F)(F)F)C#N)=S)C(NC)=O)F"]


def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('output_txt_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('start_energy', type=str, help='starting energy')
    parser.add_argument('interval', type=str, help='interval')
    parser.add_argument('experiment_for_one_molecule', type=str, help='interval')
    parser.add_argument('traj', type=str, help='traj')

    args = parser.parse_args()
    
    output_txt_file = args.output_txt_file
    start_energy = int(args.start_energy)
    interval = int(args.interval)
    experiment_per_molecule = int(args.experiment_for_one_molecule)
    traj = args.traj


    smile_list = smile_list_small
    with open(output_txt_file, 'w') as output_file:
        output_file.write("smiles,energy,trajectories,atom,realEnergy\n")
        for smile in smile_list:
            for i in range(0, experiment_per_molecule):
                current_energy = i * interval + start_energy
                output_file.write(smile + "," + str(current_energy) + "," + traj + "\n")


if __name__ == "__main__":
    main()