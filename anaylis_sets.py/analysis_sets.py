# ~/LabData/public_public/QCXMS/set1/analysis/ecom_test
#  python ~/QCxMS_Analyze/plot_box_plots.py mols_ecom_test.txt ~/QCXMS/workDir/neccessary_files/FDA-Part2/delete_precursor_20 merged_predicted_spectra.mgf ~/QCXMS/workDir/neccessary_files/ALL_GNPS_cleaned.csv ~/QCXMS/data/GNPS-SELLECKCHEM-FDA-PART2/GNPS-SELLECKCHEM-FDA-PART2.json 5 3
import os
import argparse
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import mgf
import collections
from typing import List, Tuple
import itertools
import argparse
import ast
import re
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from scipy.optimize import curve_fit
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import math


SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "mz", "intensity"]
)
def sort_array_by_array(arrayA, arrayB):
    return [x for _, x in sorted(zip(arrayB, arrayA))]

# Define the polynomial function
def func(x, a, b, c):
    return a * x**2 + b * x + c

def norm_intensity(intensity):
    return np.copy(intensity)/np.linalg.norm(intensity)

def _cosine_fast(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    allow_shift: bool,
    allow_sqrt = False,
) -> Tuple[float, List[Tuple[int, int]]]:
    precursor_charge = 1
    precursor_mass_diff = (
        #spec.precursor_mz - spec_other.precursor_mz
        0
    ) * precursor_charge
    # Only take peak shifts into account if the mass difference is relevant.
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= fragment_mz_tolerance:
        num_shifts += precursor_charge
    other_peak_index = np.zeros(num_shifts, np.uint16)
    mass_diff = np.zeros(num_shifts, np.float32)
    for charge in range(1, num_shifts):
        mass_diff[charge] = precursor_mass_diff / charge

    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - fragment_mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1
        # Match the peaks within the fragment mass window if possible.
        for cpi in range(num_shifts):
            index = 0
            other_peak_i = other_peak_index[cpi] + index
            while (
                other_peak_i < len(spec_other.mz)
                and abs(
                    peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])
                )
                <= fragment_mz_tolerance
            ):
                if not allow_sqrt:
                    peak_match_scores.append(peak_intensity * spec_other.intensity[other_peak_i])
                else:
                    peak_match_scores.append(math.sqrt(peak_intensity) * math.sqrt(spec_other.intensity[other_peak_i]))
                peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    score, peak_matches = 0.0, []
    if len(peak_match_scores) > 0:
        # Use the most prominent peak matches to compute the score (sort in
        # descending order).
        peak_match_scores_arr = np.asarray(peak_match_scores)
        peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(peak_match_idx)[peak_match_order]
        peaks_used, other_peaks_used = set(), set()
        for peak_match_score, peak_i, other_peak_i in zip(
            peak_match_scores_arr,
            peak_match_idx_arr[:, 0],
            peak_match_idx_arr[:, 1],
        ):
            if (
                peak_i not in peaks_used
                and other_peak_i not in other_peaks_used
            ):
                score += peak_match_score
                # Save the matched peaks.
                peak_matches.append((peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)

    return score, peak_matches

# return scan # <-> SpectrumTuple(mz,norm_intensity(intensity))
def get_mgf_dict(input_mgf_file, delete_precursor_switch):
    spec_dic = {}
    count = 0
    noise_reduction_flag = False
    for spectrum in mgf.read(input_mgf_file):
        #print(spectrum)
        params = spectrum.get('params')
        mz = spectrum.get('m/z array')
        intensity = spectrum.get('intensity array')
        min_intensity = min(intensity) 

        spec_dic[params['scan']] = {}

        if not delete_precursor_switch:
            spec_dic[params['scan']] = SpectrumTuple(params['pepmass'][0],mz,norm_intensity(intensity))
        else: 
            precursor = params['pepmass'][0]
            out_mz_intensity = []

            mz_out = []
            intensity_out = []
            for i in range(0, len(mz)):
                # if noise_reduction_flag and mz[i] > 50:# and intensity[i] > 3 * min_intensity:
                #     mz_out.append(mz[i])
                #     intensity_out.append(intensity[i])
                # noise_reduction_flag and mz[i] > 50:
                if mz[i] > 50:
                    mz_out.append(mz[i])
                    intensity_out.append(intensity[i])

            spec_dic[params['scan']] = SpectrumTuple(params['pepmass'][0], mz_out, norm_intensity(intensity_out))

        count = count + 1

    return spec_dic

def get_smile_peak_dict(json_smile_peak_file):
    smile_peak_dic = {}
    noise_reduction_flag = True
    with open(json_smile_peak_file, "r") as json_smile_peak_file:
        for line in json_smile_peak_file:
            if str(line) != "smiles: peaks_json, atom\n":
                smile = line.split(": ")[0]
                data = line.split(": ")[1].split(" id=")[0]
                
                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, data)

                # Convert the matches to float and group them into pairs
                coords = [[float(matches[i]), float(matches[i+1])] for i in range(0, len(matches), 2)]

                mz = [x[0] for x in coords]
                #print(mz)
                intensity = [x[1] for x in coords]
                #print(intensity)

                # remove noise:
                min_intensity = min(intensity) 
                mz_out = []
                intensity_out = []
                for i in range(0, len(mz)):
                    
                        if noise_reduction_flag and intensity[i] > 3 * min_intensity and mz[i] > 50:
                            mz_out.append(mz[i])
                            intensity_out.append(intensity[i])
                        elif not noise_reduction_flag and mz[i] > 50:
                           # print("intensity[i]", intensity[i], "3 * min_intensity", 3 * min_intensity)
                            mz_out.append(mz[i])
                            intensity_out.append(intensity[i])

                smile_peak_dic[smile] = SpectrumTuple("0", mz,norm_intensity(intensity))
    return smile_peak_dic

def get_scan_smile_dict(input_mole_file): 
    scan_smile_dict = {}
    energy_list = []
    molecule_list = []
    first_smile = ""
    with open(input_mole_file) as input_mole_file:
        scanNumber = 0
        for line in input_mole_file:
            if "smiles," not in line:
                scan_smile_dict[str(scanNumber)] = line.split(",")[0]
                if line.split(",")[0] == first_smile:
                    energy_list.append(float(line.split(",")[1]))
                    #print("same smile", line.split(",")[1])

                molecule = line.split(",")[0]
                if molecule not in molecule_list:
                    molecule_list.append(molecule)

            if scanNumber == 1:
                energy_list.append(int(line.split(",")[1]))
                first_smile = line.split(",")[0]
            scanNumber = scanNumber + 1

    return scan_smile_dict, energy_list, molecule_list

def get_smile_id_dict(json_smile_peak_file):
    smile_id_dict = {}
    #with open(input_mole_file, "r") as input_mole_file: 
    with open(json_smile_peak_file, "r") as json_smile_peak_file:
        for line in json_smile_peak_file:
            if str(line) != "smiles: peaks_json, atom\n":
                smile = line.split(": ")[0]
                #print(line)
                id = line.split(" id= ")[1].split(": ")[0]
                smile_id_dict[smile] = id

    return smile_id_dict

def get_id_energy_dict(real_energy_file):
    id_energy_dit = {}
    with open(real_energy_file, "r") as file:
        for line in file:
            if line != "	spectrum_id	collision_energy\n":
                id = line.split(",")[1]
                energy = line.split(",")[2]
                id_energy_dit[id] = energy

    return id_energy_dit
def plot_boxplot(accuracy_list, x_list, plotName="trend_boxplot.png", x_axis = "CE (eV)", y_axis = "Accuracy"):
    plt.clf()
    plt.boxplot(accuracy_list)
#    plt.title('Intensity Trend (Boxplot)')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.xticks(range(1, len(x_list) + 1), x_list)  # Set x-axis tick labels
    plt.savefig(plotName)
    plt.clf()


def plot_line_graph(accuracy_list, x_list, plotName= "trend.png"):
    plt.clf()
    plt.plot(x_list, accuracy_list)
    plt.title('intensity trend')
    plt.xlabel('collision Energy')
    plt.ylabel('accuracy')
    plt.savefig(plotName)
    plt.clf()

def plot_dual_line_graph(accuracy_list_1, accuracy_list_2, x_list, plotName="trend_dual_y.png"):
    plt.clf()

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('CE (eV)')
    ax1.set_ylabel('cosine Score')
    ax1.plot(x_list, accuracy_list_1, label='Y-axis 1 Data')
    ax1.tick_params(axis='y', labelcolor='tab:red')

    # Create a second y-axis at the bottom
    ax2 = ax1.twiny()
    ax2.set_ylabel('Y-axis 2')
    neg_accuracy_list_2 = [[element * -1 for element in sublist] for sublist in accuracy_list_2]
    ax2.plot(x_list, neg_accuracy_list_2, label='Y-axis 2 Data')
    ax2.tick_params(axis='y', labelcolor='tab:blue')

    # Set the y-axis limits for both y-axes
    ax1.set_ylim(1, 0)
    ax2.set_ylim(-1, 1)

    # Set the x-axis at y=0
    ax1.axhline(y=0, color='black', linewidth=0.5)
    ax1.annotate('x-axis', xy=(-0.025, 0.5), rotation=90, va='center', ha='right')

     # Modify tick labels for the second y-axis
    y2_ticks = [1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75, -1]
    y2_ticklabels = ['1', '0.75', '0.5', '0.25', '0', '0.25', '0.5', '0.75', '1']
    ax2.set_yticks(y2_ticks)
    ax2.set_yticklabels(y2_ticklabels)

    # Adjust subplot positioning to center the plot
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.15, right=0.85)

    plt.title('Cosine score vs. CE')
    plt.savefig(plotName)
    plt.clf()


def plot_scatter_graph(accuracy_list, x_list):
    plt.clf()
    
    for i, accuracy_values in enumerate(accuracy_list):
        plt.scatter(x_list, accuracy_values, marker='o', label=f'Set {i + 1}', s=3)
    
    plt.title('Intensity Trend')
    plt.xlabel('Collision Energy')
    plt.ylabel('Accuracy')
    #plt.legend()
    plt.savefig('trend_scatter.png')
    plt.clf()

def ana_energy_score(energy_score_list, energy_list):
    plt.clf()
    # Calculate the averages for each list
    averages = [sum(x)/len(x) for x in energy_score_list]

    # Create a bar graph of the averages
    plt.bar(energy_list, averages)

    # Add labels and title
    plt.xlabel('collision energy')
    plt.ylabel('Average')
    plt.title('Average score of Each energy')

    # Display the plot
    plt.savefig('energy_score.png')
    plt.clf()

def get_n_largest_peaks(
        smile_peak_dict,
        n_peak):

    n_largest_peaks_dict = {}
    for mol in smile_peak_dict:
        sorted_indices = sorted(range(len(smile_peak_dict[mol].intensity)), key=lambda i: smile_peak_dict[mol].intensity[i], reverse=True)
        sorted_mz = [smile_peak_dict[mol].mz[i] for i in sorted_indices[:n_peak]]
        sorted_intensity = [smile_peak_dict[mol].intensity[i] for i in sorted_indices[:n_peak]]
        n_largest_peaks_dict[mol] = SpectrumTuple(smile_peak_dict[mol].precursor_mz, sorted_mz, sorted_intensity)
    return n_largest_peaks_dict

def plot_violin(energy_score_list, energy_list):
    plt.clf()
    # Find the maximum length among the lists
    max_length = max(len(lst) for lst in energy_score_list)

    # Pad the lists with None to make them equal length
    padded_data = [lst + [None] * (max_length - len(lst)) for lst in energy_score_list]

    # Create a figure and axis
    fig, ax = plt.subplots()

    # Plot the violin plot
    sns.violinplot(data=padded_data, ax=ax)

    # Set custom x-axis tick labels
    xtick_labels = energy_list
    ax.set_xticks(range(len(xtick_labels)))
    ax.set_xticklabels(xtick_labels)

    # Set labels and title
    ax.set_xlabel('Predicting Energy')
    ax.set_ylabel('Accuracy')
    ax.set_title('Violin Plot of Energy Effect')

    # Save the plot
    plt.savefig('output/violin_plot.png')

def get_compund_info(input_smiles, input_json):
    with open(input_json, 'r') as f:
        data = json.load(f)

    # Loop through all the strings in the JSON data
    smile_dict = {}
    for item in data:        
        if 'Smiles' in str(item):  
            smile = item["Smiles"]
            Adduct = item["Adduct"]
            instrument = item["Instrument"]
            spectrum_id = item["spectrum_id"]
            USI = "mzspec:GNPS:GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE:accession:"+spectrum_id
            # if smile in smile_dict.keys():
            if smile not in smile_dict.keys() and "M+H" in Adduct:# and "Orbitrap" in instrument:     
                #print("smile ", smile_dict[smile])
                smile_dict[smile] = [USI, Adduct, instrument]

    smiles_list = []
    with open(input_smiles, 'r') as f:
        data = f.read().splitlines()
    
        for line in data:
            if "smiles" not in line:
                smiles_list.append(line.split(",")[0])

    exist_smiles_dict = {}

    for smile in smiles_list:
            exist_smiles_dict[smile]=smile_dict[smile]

    return exist_smiles_dict


def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_mole_file', type=str, help='input_mole_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('json_smile_peak_file', type=str, help='json_smile_peak_file is the converted file that has the information that each smiles in the json file with its peak and intensity')
    parser.add_argument('input_mgf_file', type=str, help='the simulation result in the mgf format.')
    parser.add_argument('real_energy_file', type=str, help=' the file stores information of experimental spetrums energy.')
    parser.add_argument('experiment_json', type=str, help=' the json file has all the information of the experimental data.')
    parser.add_argument('experiment_for_each_mol', type=str, help='is a numer that for each compound, how many experiments are performed on it')
    #parser.add_argument('tolerance', type=str, help='input_mgf_file is the mgf file from HPCC')

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    print("script_dir", script_dir)


    input_mole_file = os.path.join(script_dir, args.input_mole_file)
    #input_mole_file = args.input_mole_file
    json_smile_peak_file = os.path.join(script_dir, args.json_smile_peak_file)
    #json_smile_peak_file = args.json_smile_peak_file 
    input_mgf_file = os.path.join(script_dir, args.input_mgf_file)
    #input_mgf_file = args.input_mgf_file
    real_energy_file = os.path.join(script_dir, args.real_energy_file)
    #real_energy_file = args.real_energy_file
    experiment_json = os.path.join(script_dir, args.experiment_json)
    # experiment_json = args.experiment_json
    experiment_for_each_mol = int(args.experiment_for_each_mol)
    #tolerance = int(args.tolerance)
    n_peak = 8
    delete_precursor_switch = True

    mgf_dict = get_mgf_dict(input_mgf_file, delete_precursor_switch)
    #print(mgf_dict)

    smile_peak_dict = get_smile_peak_dict(json_smile_peak_file)
    #print("smile_peak_dict: ", smile_peak_dict)
    largest_n_peaks_spect_dict = get_n_largest_peaks(smile_peak_dict, n_peak)

    scan_smile_dict, energy_list, molecule_list = get_scan_smile_dict(input_mole_file)
    #print("energy_list: ", energy_list)

    id_energy_dict = get_id_energy_dict(real_energy_file)
    #print("id_energy_dit: ", id_energy_dit)

    smile_id_dict = get_smile_id_dict(json_smile_peak_file)
    #print("smile_id_dict: ", smile_id_dict)

    smile_info_dict = get_compund_info(input_mole_file, experiment_json)

    # structure
    # mgf_dict:                             scan# <=> mgf's peak_intensity
    # scan_smile_dict: scan# <=> smile        |
    # smile_peak_dict:                      smile <=> json's peak_intensity
    # smile_id_dict:    smile <=> id           |
    # id_energy_dit:                         id   <=> energy
   
    score_list=[]
    explained_list = []
    explained_num = []
    energy_score_list = []
    sqrt_energy_score_list = []
    row_matched_intensity_list = []
    row_matched_num_list = []
    simulation_num_list = []
    simulation_explained_list = []
    table_list = []
    for i in range(0, experiment_for_each_mol):
        energy_score_list.append([])
        explained_list.append([])
        explained_num.append([])
        sqrt_energy_score_list.append([])
        row_matched_intensity_list.append([])
        row_matched_num_list.append([])
        simulation_num_list.append([])
        simulation_explained_list.append([])

    missing_counter = 0
    smile_groundTruth_energy_dict = {}
    each_mol_matched_list = []
    simulated_mol_matched_list = []
    matched_percentage_pre = 0
    simulated_matched_percentage = 0
    # scan starts from 1:
    for i in range(1, len(scan_smile_dict)+1) :
        if i % len(energy_list) == 1:
            each_mol_matched_list = []
            simulated_mol_matched_list = []
        
        row_matched_num = 0
        experimental_peaks_num = 0
        simulated_peaks_num = 0
        explained_intensity = 0
        smile = scan_smile_dict[str(i)]
        if str(i) in mgf_dict.keys():        
            counter = 1
            
            if smile+"_" + str(counter) not in smile_peak_dict:
                continue
            if smile+"_" + str(counter) in smile_groundTruth_energy_dict.keys():                   
                smile_groundTruth_energy_dict[smile] = id_energy_dict[smile_id_dict[smile]]
            else:
                smile_groundTruth_energy_dict[smile] = id_energy_dict[smile_id_dict[smile+"_" + str(counter)]]

            spec1=mgf_dict[str(i)]
            spec2=smile_peak_dict[smile+"_" + str(counter)]
            score, matched_peaks=_cosine_fast(spec1,spec2,0.1,False)
            spec1_sqrt = SpectrumTuple(spec1.precursor_mz,spec1.mz,norm_intensity([math.sqrt(i) for i in spec1.intensity]))
            spec2_sqrt = SpectrumTuple(spec2.precursor_mz,spec2.mz,norm_intensity([math.sqrt(i) for i in spec2.intensity]))
            score_sqrt, matched_peaks_sqrt=_cosine_fast(spec1_sqrt,spec2_sqrt,0.1,False)
            for tuple in matched_peaks:
                if spec2.intensity[tuple[1]] not in each_mol_matched_list:
                    each_mol_matched_list.append(spec2.intensity[tuple[1]])
                    simulated_mol_matched_list.append(spec1.intensity[tuple[0]])

            raw_matched_intensity = 0
            simulated_matched_intensity = 0
            for tuple in matched_peaks:
                raw_matched_intensity +=  spec2.intensity[tuple[1]]
                simulated_matched_intensity +=  spec1.intensity[tuple[0]]

            matched_intensity = sum([j for j in each_mol_matched_list])
            matched_num = len(each_mol_matched_list)
            all_intensity = sum(spec2.intensity)
            matched_percentage = matched_intensity / all_intensity * 100
            raw_matched_percentage = raw_matched_intensity / sum(spec2.intensity) * 100
            simulated_matched_percentage = simulated_matched_intensity / sum(spec1.intensity) * 100
            matched_percentage_pre = matched_percentage
            simulated_matched_percentage_pre = simulated_matched_percentage * 100

            if len(matched_peaks) <= 3:
                score = 0

            row_matched_num = len(matched_peaks) 
            experimental_peaks_num = len(spec2.intensity)
            explained_intensity = matched_percentage
            simulated_peaks_num = len(spec1.intensity)

            energy_score_list[i % experiment_for_each_mol].append(score * 100)
            sqrt_energy_score_list[i % experiment_for_each_mol].append(score_sqrt)
            row_matched_intensity_list[i % experiment_for_each_mol].append(raw_matched_percentage)
            row_matched_num_list[i % experiment_for_each_mol].append(row_matched_num/ len(spec2.intensity)* 100)
            explained_list[i % experiment_for_each_mol].append(matched_percentage)
            explained_num[i % experiment_for_each_mol].append(matched_num/len(spec2.intensity) * 100)
            simulation_num_list[i % experiment_for_each_mol].append(len(matched_peaks) / len(spec1.intensity) * 100)
            simulation_explained_list[i % experiment_for_each_mol].append(simulated_matched_percentage)

            score_list.append(score)

            
        else:
            print("do not have this scan: ", i)
            score = 0
            score_list.append(score)
            energy_score_list[i % experiment_for_each_mol].append(score)
            sqrt_energy_score_list[i % experiment_for_each_mol].append(score_sqrt)
            row_matched_intensity_list[i % experiment_for_each_mol].append(0)
            row_matched_num_list[i % experiment_for_each_mol].append(0)
            explained_list[i % experiment_for_each_mol].append(0)
            explained_num[i % experiment_for_each_mol].append(0)
            simulation_num_list[i % experiment_for_each_mol].append(0)
            simulation_explained_list[i % experiment_for_each_mol].append(0)
            
            missing_counter = missing_counter + 1
            if smile+"_" + str(counter) in smile_groundTruth_energy_dict.keys():                   
                smile_groundTruth_energy_dict[smile] = id_energy_dict[smile_id_dict[smile]]
            else:
                smile_groundTruth_energy_dict[smile] = id_energy_dict[smile_id_dict[smile+"_" + str(counter)]]

        #items scan, smile, energy, num_of_peaks_of_experimental_spectrum, number_of_peaks_in_simulated_spectrum, explained_intensity
        item = [i, scan_smile_dict[str(i)], energy_list[(i-1)%len(energy_list)], experimental_peaks_num, simulated_peaks_num, explained_intensity]
        table_item = item[:2] + smile_info_dict[scan_smile_dict[str(i)]] + item[2:]
        table_list.append(table_item)
           

    ave_list = []
    sorted_list = []
    sqrt_sorted_list = []
    explained_intensity_list = []
    raw_explained_intensity_list = []
    raw_explained_num_list = []
    explained_num_list = []
    simulated_raw_explained_num_list = []
    simulation_explained_intensity_list = []
    for i in range(1, experiment_for_each_mol + 1):
        if i == experiment_for_each_mol:
            ave_list.append(sum(energy_score_list[0]) / len(energy_score_list[0]))
            sorted_list.append((energy_score_list[0]))
            sqrt_sorted_list.append((sqrt_energy_score_list[0]))
            raw_explained_intensity_list.append(row_matched_intensity_list[0])
            raw_explained_num_list.append(row_matched_num_list[0])
            explained_intensity_list.append((explained_list[0]))
            explained_num_list.append((explained_num[0]))
            simulated_raw_explained_num_list.append(simulation_num_list[0])
            simulation_explained_intensity_list.append((simulation_explained_list[0]))

        else: 
            ave_list.append(sum(energy_score_list[i]) / len(energy_score_list[i]))
            sorted_list.append((energy_score_list[i]))
            sqrt_sorted_list.append((sqrt_energy_score_list[i]))
            raw_explained_intensity_list.append(row_matched_intensity_list[i])
            raw_explained_num_list.append(row_matched_num_list[i])
            explained_intensity_list.append((explained_list[i]))
            explained_num_list.append((explained_num[i]))
            simulated_raw_explained_num_list.append(simulation_num_list[i])
            simulation_explained_intensity_list.append(simulation_explained_list[i])


    sorted_list_T = np.array(sorted_list).T
    raw_explained_intensity_list_T = np.array(raw_explained_intensity_list).T
    sorted_raw_explained_intensity_list_T = [sorted(sublist) for sublist in raw_explained_intensity_list_T]
    sorted_raw_explained_intensity_list = np.array(sorted_raw_explained_intensity_list_T).T
    raw_explained_num_list_T = np.array(raw_explained_num_list).T
    sorted_raw_explained_num_list_T = [sorted(sublist) for sublist in raw_explained_num_list_T]
    sorted_raw_explained_num_list = np.array(sorted_raw_explained_num_list_T).T
    best_cosine_for_each_molecule_index = np.argmax(sorted_list_T, axis=1)
    best_intensity__each_molecule_index = np.argmax(raw_explained_intensity_list_T, axis=1)
    #best_energy_values = [[energy_list[index]] for index in best_cosine_for_each_molecule_index]
    best_energy_values = [energy_list[index] for index in best_cosine_for_each_molecule_index]
    best_intensity__values = [energy_list[index] for index in best_intensity__each_molecule_index]

    temp_list = [i for i in range(1, len(energy_list) + 1)]
    #plot_dual_line_graph(sorted_raw_explained_intensity_list, sorted_raw_explained_num_list, temp_list, "predicted_numbers.png")
    #plot_line_graph(sorted_list, energy_list, "cosine_line.png")
    #plot_line_graph(sqrt_sorted_list, energy_list, "sqrt_cosine_line.png")
    plot_boxplot(sorted_list, energy_list, "output/cosine_score.png", "CE (eV)", "Cosine Score")
    plot_boxplot(raw_explained_intensity_list, energy_list, "output/percentage_intensity_boxplot.png", "CE (eV)", "Explained Intensity")
    plot_boxplot(raw_explained_num_list, energy_list, "output/percentage_num_boxplot.png", "CE (eV)", "% of explained peaks")
    molecule_list = [i for i in range(0, 1)]
    plot_boxplot(best_energy_values, molecule_list, "output/each_mole_boxplot.png", "Molecule", "CE (eV)")
    plot_boxplot(best_intensity__values, molecule_list, "output/each_mole_boxplot_intensity.png", "molecule", "CE (eV)")
    plot_boxplot(simulated_raw_explained_num_list, energy_list, "output/simulated_percentage_num_boxplot.png", "CE (eV)", "% of explained peaks")
    plot_boxplot(simulation_explained_intensity_list, energy_list, "output/simulated_percentage_intensity_boxplot.png", "CE (eV)", "Explained Intensity")

    # write the table
    with open("output/table.txt", "w") as table_output:
        table_output.write("scan,smile,USI,Adduct,instrument,energy,num_of_peaks_of_experimental_spectrum,number_of_peaks_in_simulated_spectrum, explained_intensity\n")
        for item in table_list:
            table_output.write(','.join(map(str, item)) + "\n")

    #plot_scatter_graph(np.array(explained_intensity_list).T, energy_list)
    
    # --------------------------below has not been cleaned -----------------------------------
    #ana_energy_score(energy_score_list, energy_list)
    plot_violin(energy_score_list, energy_list)
    plt.clf()

    sns.heatmap(sorted_list_T, cmap='YlGnBu', annot=False, xticklabels=energy_list)
    plt.title('Accuracy of Collision Energy')
    plt.xlabel('collision Energy')
    plt.ylabel('molecules')
    plt.savefig('output/spec.png')


    #analysis n_largest_intensity match 
    num_of_matches = []
    for i in range(1, len(scan_smile_dict)+1) :
        #print(mgf_dict.keys())
        smile = scan_smile_dict[str(i)]
        if str(i) in mgf_dict.keys():        
            counter = 1
            if smile+"_" + str(counter) not in smile_peak_dict:
                continue
            spec1=mgf_dict[str(i)]
            spec2=largest_n_peaks_spect_dict[smile+"_" + str(counter)]
            score, matched_peaks=_cosine_fast(spec1,spec2,0.1,False)
            num_of_matches.append(matched_peaks)

            n_peak_over_1 = len(matched_peaks) / n_peak * 100
            #print(f"matching peaks out of {n_peak}: {len(matched_peaks)}, maching percenation is: {n_peak_over_1}%")

    molecule_energy_score_list = np.array(sorted_list).T
    #print(molecule_energy_score_list)
    max_molecule_energy_score_list = []
    max_score_list = []
    useful_molecule_index = []
    #print("max_molecule_energy_score_list: ", max_molecule_energy_score_list)
    counter = 0
    for i in molecule_energy_score_list:
        max_score_list.append(max(i))
        max_molecule_energy_score_list.append(energy_list[np.argmax(i)])
        #store the molecule that has the highest score > shreshold
        if max(i) > 0.7:
            useful_molecule_index.append(counter)
            #print("useful molecule: ", molecule_list[counter])
        counter += 1

    # not a function, but to produce the not good predictions
    counter = 0
    for i in molecule_energy_score_list:
        if max(i) < 0.7:
            print("bad molecule: ", molecule_list[counter])
        counter += 1
  
    print("useful_molecule_index: ", useful_molecule_index)
    print("max_molecule_energy_score_list", max_molecule_energy_score_list)
    print("smile_groundTruth_energy_dict: ", smile_groundTruth_energy_dict)
    print("max_score_list: ", max_score_list)
    for i in max_score_list:
        print(i)
    groundTruth_energy_list = [float(x) for x in smile_groundTruth_energy_dict.values()]

    max_molecule_energy_score_list = np.array(max_molecule_energy_score_list).reshape(-1, 1)
    poly = PolynomialFeatures(degree=1)
    x_poly = poly.fit_transform(np.array(max_molecule_energy_score_list).reshape(-1, 1))
    regressor = LinearRegression()
    regressor.fit(x_poly, np.array(groundTruth_energy_list))

    y_pred = regressor.predict(x_poly)

    # Calculate the R-squared value
    r2 = r2_score(np.array(groundTruth_energy_list), y_pred)
    print("R-squared value:", r2)

    # Print the regression results
    print("Coefficients: ", regressor.coef_)
    print("Intercept: ", regressor.intercept_)
    
    #only use those useful useful_molecule to do regression
    useful_molecule_energy_score_list = []
    useful_groundTruth_energy_list = []
    for i in useful_molecule_index:
        print("useful_molecule_index: ", i)
        useful_molecule_energy_score_list.append(max_molecule_energy_score_list[i])
        useful_groundTruth_energy_list.append(groundTruth_energy_list[i])

    print("useful_molecule_energy_score_list", useful_molecule_energy_score_list)
    print("useful_groundTruth_energy_list: ", useful_groundTruth_energy_list)
    
    useful_molecule_energy_score_list = np.array(useful_molecule_energy_score_list).reshape(-1, 1)
    poly = PolynomialFeatures(degree=1)
    x_poly = poly.fit_transform(np.array(useful_molecule_energy_score_list).reshape(-1, 1))
    regressor = LinearRegression()
    regressor.fit(x_poly, np.array(useful_groundTruth_energy_list))

    y_pred = regressor.predict(x_poly)



    # Calculate the R-squared value
    r2 = r2_score(np.array(useful_groundTruth_energy_list), y_pred)
    print("useful R-squared value:", r2)

    # Print the regression results
    print("useful Coefficients: ", regressor.coef_)
    print("useful Intercept: ", regressor.intercept_)


    

# add legend
plt.legend()



if __name__ == "__main__":
    main()
