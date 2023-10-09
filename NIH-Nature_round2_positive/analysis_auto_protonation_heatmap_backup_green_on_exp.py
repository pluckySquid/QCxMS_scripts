# ~/LabData/public_public/QCXMS/GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE/result/first_25_atoms_result$
# python ~/QCXMS/scripts/NIH-Nature_round2_positive/analysis_auto_protonation_heatmap_backup_green_on_exp.py mols_first_25_atoms.txt ~/QCXMS/workDir/neccessary_files/NIH-Nature_round2_positive/NIH_Nature_round2_positive_delete_precursor_17 merged_predicted_spectra.mgf ~/QCXMS/workDir/neccessary_files/ALL_GNPS_cleaned.csv 1 3 > debug.txt
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
from scipy.optimize import minimize
import csv


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

    #print("len other_peak_index: ", len(other_peak_index))
    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        #print("peak_index: ", peak_index)
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - fragment_mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1

        #print("other_peak_index: ", other_peak_index)
        #print("num_shifts: ", num_shifts)
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
                peak_match_scores.append(
                    peak_intensity * spec_other.intensity[other_peak_i]
                )
                peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    #print("peak_match_scores: ", peak_match_scores)
    #print("peak_match_idx: ", peak_match_idx)
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
                #print("(peak_i, other_peak_i): ", (peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)
            #else:
                #print("used peaks ", (peak_i, other_peak_i))

    return score, peak_matches

# Define the cosine similarity score function
def cosine_score_function(weights, lists):
    weights = np.array(weights, dtype=float)
    # print("weights, ", weights)
    combined_list = np.sum([w * np.array(lst) for w, lst in zip(weights, lists[:-1])], axis=0)
    # print("combined_list: ", combined_list)
    # print("experiment list: ", combined_list[-1])

    predict_mz = [x[0] for x in combined_list]
    predict_intensity = [x[1] for x in combined_list]
    predict_spec = SpectrumTuple("0", predict_mz, norm_intensity(predict_intensity))
                                 
    experiment_mz = [x[0] for x in  lists[-1]]
    experiment_intensity = [x[1] for x in  lists[-1]]
    experiment_spec = SpectrumTuple("0", experiment_mz, norm_intensity(experiment_intensity))                             
    
    # print("experiment_spec: ", experiment_spec)
    cosine_score = cal_cosine_score(predict_spec, experiment_spec, 0.02, True)
    # print("cosine_score: ", cosine_score)

    cosine_fast_score = _cosine_fast(predict_spec, experiment_spec, 0.02, True)
    #print("cosine_fast_score: ", cosine_fast_score)

    return -cosine_score  # Negate the score for minimization

def cal_cosine_score(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    allow_shift: bool,
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

    #print("len other_peak_index: ", len(other_peak_index))
    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        #print("peak_index: ", peak_index)
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - fragment_mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1

        #print("other_peak_index: ", other_peak_index)
        #print("num_shifts: ", num_shifts)
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
                peak_match_scores.append(
                    peak_intensity * spec_other.intensity[other_peak_i]
                )
                peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    #print("peak_match_scores: ", peak_match_scores)
    #print("peak_match_idx: ", peak_match_idx)
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
                #print("(peak_i, other_peak_i): ", (peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)
            #else:
                #print("used peaks ", (peak_i, other_peak_i))

    return score

# return scan # <-> SpectrumTuple(mz,norm_intensity(intensity))
def get_mgf_dict(input_mgf_file, delete_precursor_switch):
    spec_dic = {}
    count = 0
    noise_reduction_flag = False
    for spectrum in mgf.read(input_mgf_file):
        #print(spectrum)
        params = spectrum.get('params')
        mz = spectrum.get('m/z array')
        intensity = [math.sqrt(x) for x in spectrum.get('intensity array')]
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
                if mz[i] > 50 and intensity[i] > 0.1 and mz[i] < params['pepmass'][0] - 17:
                    mz_out.append(mz[i])
                    intensity_out.append(intensity[i])

            #print("np.array(mz_out): ", mz_out)
            spec_dic[params['scan']] = SpectrumTuple(params['pepmass'][0], mz_out, norm_intensity(intensity_out))
        #print(spec_dic[params['scan']])
        count = count + 1

    return spec_dic

def get_smile_peak_dict(json_smile_peak_file):
    smile_peak_dic = {}
    #with open(input_json_smile_file, "r") as input_json_smile_file: 
    with open(json_smile_peak_file, "r") as json_smile_peak_file:
        for line in json_smile_peak_file:
            if str(line) != "smiles: peaks_json, atom\n":
                smile = line.split(": ")[0]
                print("smile:", smile)
                data = line.split(": ")[1].split(" id=")[0]
                
                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, data)

                # Convert the matches to float and group them into pairs
                coords = [[float(matches[i]), float(matches[i+1])] for i in range(0, len(matches), 2)]

                mz = [x[0] for x in coords]
                #print(mz)
                intensity = [x[1] for x in coords]
                normlised_intensity = norm_intensity(intensity)
                print("intensity: ", normlised_intensity)
                denoised_intensity = [math.sqrt(x) for x in normlised_intensity if x > 0.02]
                indexes = [index for index, value in enumerate(normlised_intensity) if value > 0.02]
                mz_indexes = [mz[index] for index in indexes]

                #print("denoised_intensity: ", denoised_intensity)

                smile_peak_dic[smile] = SpectrumTuple("0", mz_indexes,denoised_intensity)
    return smile_peak_dic

def get_scan_smile_dict(input_json_smile_file): 
    scan_smile_dict = {}
    energy_list = []
    molecule_list = []
    first_smile = ""
    with open(input_json_smile_file) as input_json_smile_file:
        scanNumber = 0
        for line in input_json_smile_file:
            if "smiles," not in line:
                scan_smile_dict[str(scanNumber)] = line.split(",")[0]
                if line.split(",")[0] == first_smile:
                    energy_list.append(int(line.split(",")[1]))
                    # print("same smile", line.split(",")[1])

                molecule = line.split(",")[0]
                if molecule not in molecule_list:
                    molecule_list.append(molecule)

            if scanNumber == 1:
                energy_list.append(float(line.split(",")[1]))
                first_smile = line.split(",")[0]
            scanNumber = scanNumber + 1

    # formanual protonation
    #energy_list = [70]
    # print(energy_list)

    return scan_smile_dict, energy_list, molecule_list

def get_smile_id_dict(json_smile_peak_file):
    smile_id_dict = {}
    #with open(input_json_smile_file, "r") as input_json_smile_file: 
    with open(json_smile_peak_file, "r") as json_smile_peak_file:
        for line in json_smile_peak_file:
            if str(line) != "smiles: peaks_json, atom\n":
                smile = line.split(": ")[0]
                #print(line)
                id = line.split(" id= ")[1].split(": ")[0]
                smile_id_dict[smile] = id
    #print(smile_peak_dic)
    return smile_id_dict

def get_id_energy_dict(real_energy_file):
    id_energy_dit = {}
    with open(real_energy_file, "r") as file:
        for line in file:
            if line != "	spectrum_id	collision_energy\n":
                id = line.split(",")[1]
                energy = line.split(",")[2]
                id_energy_dit[id] = energy
                #print(id, energy)
    return id_energy_dit

def plot_line_graph(accuracy_list, x_list):
    plt.clf()
    # print(accuracy_list)
    plt.plot(x_list, accuracy_list)
    plt.title('Accuracy trend')
    plt.xlabel('collision Energy')
    plt.ylabel('accuracy')
    plt.savefig('trend.png')
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
        sorted_intensity = [math.sqrt(x) for x in sorted_intensity]

        peak_match_scores_arr = np.asarray(sorted_mz)
        peak_match_order = np.argsort(peak_match_scores_arr)
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(sorted_intensity)[peak_match_order]

        sorted_mz = peak_match_scores_arr
        sorted_intensity = peak_match_idx_arr
        
        n_largest_peaks_dict[mol] = SpectrumTuple(smile_peak_dict[mol].precursor_mz, sorted_mz, sorted_intensity)
        # print("n_largest_peaks_dict[mol]:", mol, n_largest_peaks_dict[mol])
    return n_largest_peaks_dict

def rank_by_intensity(
        smile_peak_dict,
        n_peak):

    n_largest_peaks_dict = {}
    for mol in smile_peak_dict:
        sorted_indices = sorted(range(len(smile_peak_dict[mol].intensity)), key=lambda i: smile_peak_dict[mol].intensity[i], reverse=True)
        sorted_mz = [smile_peak_dict[mol].mz[i] for i in sorted_indices[:n_peak]]
        sorted_intensity = [smile_peak_dict[mol].intensity[i] for i in sorted_indices[:n_peak]]
        sorted_intensity = [math.sqrt(x) for x in sorted_intensity]

        n_largest_peaks_dict[mol] = SpectrumTuple(smile_peak_dict[mol].precursor_mz, sorted_mz, sorted_intensity)
        print("rank_by_intensity:", mol, n_largest_peaks_dict[mol])
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
    plt.savefig('violin_plot.png')

def plot_boxplot(accuracy_list, x_list, plotName="trend_boxplot.png", x_axis = "CE (eV)", y_axis = "Accuracy"):
    plt.clf()
    plt.boxplot(accuracy_list)
#    plt.title('Intensity Trend (Boxplot)')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.xticks(range(1, len(x_list) + 1), x_list)  # Set x-axis tick labels
    plt.savefig(plotName)
    plt.clf()

def plot_barplot(describing_dict, experimental_smile_intensity_dict, plotName="explain_proteomer_barplot.png"):
# [experimental_smile_intensity_dict[smile][1][j] for j in second_elements]
# our goal is to see if a peak is explained by one or more proteomers
    #print("describing_dict", describing_dict)

    mol_single_proteomer_expalin_intensity_dict = {}
    mol_multi_proteomer_expalin_intensity_dict = {}
    mol_all_proteomer_expalin_intensity_dict = {}
    mol_cannot_explain_proteomer_expalin_intensity_dict = {}
    mol_total_intensity_dict = {}
    #interate over each smile in the experimental data
    for experiment_smile in describing_dict:    
        # check if there is only one proteomer, if so, this smile cannot count as a multi-proteomer
        #if len(describing_dict[experiment_smile])  == 1:
        #    continue   
        mol_single_proteomer_expalin_intensity_dict[experiment_smile] = []
        mol_multi_proteomer_expalin_intensity_dict[experiment_smile] = []
        mol_all_proteomer_expalin_intensity_dict[experiment_smile] = []
        mol_cannot_explain_proteomer_expalin_intensity_dict[experiment_smile] = []
        # #interate over each peaks in the experimental data:
        print("experimental_smile_intensity_dict[experiment_smile]", experimental_smile_intensity_dict[experiment_smile])
        for peak_index in range(0, len(experimental_smile_intensity_dict[experiment_smile][1])):
            proteomer_expalin_intensity = 0
            exist_num_counter = 0
            #interate over each proteomer in the experimental data
            for proteomer_index in range(0, len(describing_dict[experiment_smile])):           
                experiment_exist_peaks_elements = [t[1] for t in describing_dict[experiment_smile][proteomer_index]]
                # check if the experiment peak is in this proteomer simulation:
                print("experiment_exist_peaks_elements", experiment_exist_peaks_elements, "peak_index", peak_index)
                if peak_index in experiment_exist_peaks_elements:
                    # this peak has been predicted once (again)
                    exist_num_counter += 1
            
            print("exist_num_counter", exist_num_counter, "len(describing_dict[experiment_smile])", len(describing_dict[experiment_smile]))
            if exist_num_counter == 0:
                mol_cannot_explain_proteomer_expalin_intensity_dict[experiment_smile].append(experimental_smile_intensity_dict[experiment_smile][1][peak_index])
            elif exist_num_counter == 1:
                mol_single_proteomer_expalin_intensity_dict[experiment_smile].append(experimental_smile_intensity_dict[experiment_smile][1][peak_index])
            elif exist_num_counter == len(describing_dict[experiment_smile]):
                mol_all_proteomer_expalin_intensity_dict[experiment_smile].append(experimental_smile_intensity_dict[experiment_smile][1][peak_index])
            else:
                mol_multi_proteomer_expalin_intensity_dict[experiment_smile].append(experimental_smile_intensity_dict[experiment_smile][1][peak_index])
            # proteomer_expalin_intensity += sum(experimental_smile_intensity_dict[experiment_smile][1][i] for i in experiment_exist_peaks_elements)
            # if experiment_smile not in mol_expalin_intensity_dict:
            #     mol_expalin_intensity_dict[experiment_smile] = proteomer_expalin_intensity
            # else:
            #     mol_expalin_intensity_dict[experiment_smile].append(proteomer_expalin_intensity)
        
        total_intensity = sum (experimental_smile_intensity_dict[experiment_smile][1])
        if experiment_smile not in mol_total_intensity_dict:
                mol_total_intensity_dict[experiment_smile] = total_intensity
        
    print("mol_single_proteomer_expalin_intensity_dict", mol_single_proteomer_expalin_intensity_dict)      
    print("mol_multi_proteomer_expalin_intensity_dict", mol_multi_proteomer_expalin_intensity_dict)      
    print("mol_all_proteomer_expalin_intensity_dict", mol_all_proteomer_expalin_intensity_dict)   
    print("mol_cannot_explain_proteomer_expalin_intensity_dict", mol_cannot_explain_proteomer_expalin_intensity_dict)
    # make the array for plot:
    single_bar = []
    multi_bar = []
    all_bar = []
    cant_explain_bar = []
    overall_total_intensity = []
    explained_intensity_list = []
    for experiment_smile in describing_dict:
        if experiment_smile in mol_single_proteomer_expalin_intensity_dict:
            print("experiment_smile in the mol_single_proteomer_expalin_intensity_dict", experiment_smile)
            single_bar.append(sum(mol_single_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100)
            multi_bar.append(sum(mol_multi_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100)
            all_bar.append(sum(mol_all_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100)
            cant_explain_bar.append(sum(mol_cannot_explain_proteomer_expalin_intensity_dict[experiment_smile])/ mol_total_intensity_dict[experiment_smile] * 100)
            overall_total_intensity.append(mol_total_intensity_dict[experiment_smile])
            explained_intensity_list.append(sum(mol_single_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100 + sum(mol_multi_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100 + sum(mol_all_proteomer_expalin_intensity_dict[experiment_smile])/mol_total_intensity_dict[experiment_smile] * 100)
        else:
            single_bar.append(0)
            multi_bar.append(0)
            all_bar.append(0)
            cant_explain_bar.append(1)
            overall_total_intensity.append(0)
            explained_intensity_list.append(0)
       # mol_total_intensity_dict[experiment_smile])


    

    x_axis = range(0, len(describing_dict))   

    # Create the stacked bar plot
    
    plt.bar(x_axis, single_bar, color='r', label='Single Proteomer')
    plt.bar(x_axis, multi_bar, bottom=single_bar, color='b', label='Multi Proteomer')
    plt.bar(x_axis, all_bar, bottom=[single + multi for single, multi in zip(single_bar, multi_bar)], color='g', label='All Proteomer')
    plt.bar(x_axis, cant_explain_bar, bottom=[single + multi + all_ for single, multi, all_ in zip(single_bar, multi_bar, all_bar)], color='gray', label='Missed Explained intensities')


    # Add legend
    plt.legend()

    # Add labels, title, etc.
    plt.xlabel('molecules ID')
    plt.ylabel('explained intensities (%)')
    plt.title('explained intensities by molecules')

    plt.tight_layout()
    plt.savefig(plotName)




    # calculate the overall
    plt.clf()
    overall_single_bar = []
    overall_multi_bar = []
    overall_all_bar = []
    overall_cant_explain_bar = []
    total_explained_intensity = []
    for experiment_smile in describing_dict:
        overall_single_bar.append(sum(single_bar)/len(describing_dict))
        overall_multi_bar.append(sum(multi_bar)/ len(describing_dict))
        overall_all_bar.append(sum(all_bar) /len(describing_dict))
        total_explained_intensity.append(sum(single_bar)/len(describing_dict) + sum(multi_bar)/ len(describing_dict) + sum(all_bar) /len(describing_dict))
        overall_cant_explain_bar.append(sum(cant_explain_bar)/ len(describing_dict))
    
    x_axis = range(0, 1)

    plt.bar(x_axis, overall_single_bar, color='r', label='Single Proteomer')
    plt.bar(x_axis, overall_multi_bar, bottom=overall_single_bar, color='b', label='Multi Proteomer')
    plt.bar(x_axis, overall_all_bar, bottom=[single + multi for single, multi in zip(overall_single_bar, overall_multi_bar)], color='g', label='All Proteomer')
    plt.bar(x_axis, overall_cant_explain_bar, bottom=[single + multi + all_ for single, multi, all_ in zip(overall_single_bar, overall_multi_bar, overall_all_bar)], color='gray', label='Missed Explained intensities')

    # Add labels, title, etc.
    plt.xlabel('Overall Molecule')
    plt.ylabel('Overall explained intensities (%)')
    plt.title('Overall explained intensities by molecules')

    plt.tight_layout()
    plt.savefig("Overall " +plotName)


    #plot box plot
    print("explained_intensity_list", explained_intensity_list)
    plot_boxplot(explained_intensity_list, x_axis, plotName="explained_intensity_box_plot", x_axis = "CE (eV)", y_axis = "Accuracy")
       


def plot_discription_rate(figure, matching_peak, energy_list, line, color, label, location):
    # print("matching_peak: ", matching_peak)
    # print("energy_list: ", energy_list)
    figure.plot(energy_list, matching_peak, line, color=color, label=label)
    figure.set_title('discrption rate')
    figure.set_xlabel('molecule')
    figure.set_ylabel('discrption rate')
    figure.legend(loc=location)
    #plt.clf()

def plot_discription_function(describing_dict, predicted_mgf_dict, smile_peak_dict):
    result_dict = {}
    each_smile_dict = {}
    experimental_smile_intensity_dict = {}
    predicted_smile_intensity_dict = {}
    
    for experiment_smile in describing_dict:
        proteomer_counter = 1
        for proteomer in describing_dict[experiment_smile]:
            result_dict[experiment_smile + "_"  + str(proteomer_counter)] = []
            # print("proteomer ", proteomer)
            # print("experiment_smile ", smile_peak_dict[experiment_smile+"_1"].mz)
            predicted_index = []
            for i in proteomer:
                predicted_index.append(i[1])
            # print("predicted_index: ", predicted_index)
            for i in range(0, len(smile_peak_dict[experiment_smile+"_1"].mz)):
                if i in predicted_index:
                    result_dict[experiment_smile + "_"  + str(proteomer_counter)].append(1)
                else:
                    result_dict[experiment_smile + "_"  + str(proteomer_counter)].append(0)
            if experiment_smile not in each_smile_dict.keys():
                each_smile_dict[experiment_smile] = [result_dict[experiment_smile + "_"  + str(proteomer_counter)]]
                experimental_smile_intensity_dict[experiment_smile] = (smile_peak_dict[experiment_smile+"_1"].mz, smile_peak_dict[experiment_smile+"_1"].intensity)
                predicted_smile_intensity_dict[experiment_smile] = [(predicted_mgf_dict[experiment_smile][proteomer_counter-1].mz, predicted_mgf_dict[experiment_smile][proteomer_counter-1].intensity)]
                #proteomer_intensity_dict[experiment_smile] = [predicted_mgf_dict[experiment_smile][proteomer_counter-1].intensity]
            else:
                each_smile_dict[experiment_smile].append(result_dict[experiment_smile + "_"  + str(proteomer_counter)])
                #experimental_smile_intensity_dict[experiment_smile].append((smile_peak_dict[experiment_smile+"_1"].mz, smile_peak_dict[experiment_smile+"_1"].intensity))
                predicted_smile_intensity_dict[experiment_smile].append((predicted_mgf_dict[experiment_smile][proteomer_counter-1].mz, predicted_mgf_dict[experiment_smile][proteomer_counter-1].intensity))
                #proteomer_intensity_dict[experiment_smile].append([predicted_mgf_dict[experiment_smile][proteomer_counter-1].intensity])

            proteomer_counter += 1

    #print("result_dict: ", proteomer_intensity_dict)
    # result_dict = {
    #     'list1': [0, 1, 0, 1],
    #     'list2': [1, 0, 0, 1, 0],
    #     'list3': [1, 1, 1],
    #     'list4': [0, 0, 1, 1, 0],
    # }
    # Convert dictionary values to a 2D numpy array
    print("result_dict", result_dict)
    matrix = np.zeros((len(result_dict), max(map(len, result_dict.values()))))

    for i, (key, value) in enumerate(result_dict.items()):
        matrix[i, :len(value)] = value
    # print("matrix: ", matrix)

    # Create heatmap
    plt.imshow(matrix, cmap='Blues', aspect='auto', interpolation='nearest')

    # Set tick labels
    plt.xticks(range(matrix.shape[1]), range(1, matrix.shape[1] + 1))
    plt.yticks(range(matrix.shape[0]), result_dict.keys())

    # Show colorbar
    plt.colorbar()

    # Save the figure
    plt.savefig('heatmap.png')
    plt.clf()

    #for proteomer_counter in result_dict plot a png:
    for smile in each_smile_dict:
        # print("each_smile_dict[smile]: ", each_smile_dict[smile])
        matrix = np.zeros((len(each_smile_dict[smile]) + 1, max(map(len, each_smile_dict[smile]))))

        for i, value in enumerate(each_smile_dict[smile]):
            matrix[i, :len(value)] = value
            # print("i: ", i)
            # print("value:", value)
        matrix[i+1, :len(each_smile_dict[smile][0])] = [1] * len(each_smile_dict[smile][0])
        # Create heatmap
        plt.imshow(matrix, cmap='Blues', aspect='auto', interpolation='nearest')

        # Set tick labels
        y_ticks = []
        for m in range(0, len(each_smile_dict[smile])):
            y_ticks.append(str(smile) +"_"  + str(m+1))
        y_ticks.append("experimental")
        plt.xticks(range(matrix.shape[1]), range(1, matrix.shape[1] + 1))
        plt.yticks(range(matrix.shape[0]), y_ticks)
        #debug
        # print("matrix.shape[0]", matrix.shape[0])
        # print("matrix.shape[1]", matrix.shape[1])

        # Show colorbar
        plt.colorbar()

        # Save the figure
        try:
            plt.savefig(str(smile)+'.png')
        except:
            print("weird \\ problem")
        plt.clf()

    proteomer_intensity_dict = {}
    #To plot a mgf plot:
    for smile in predicted_smile_intensity_dict:
        smile_list = predicted_smile_intensity_dict[smile]
        #print("smile: ", smile)
        #print("smile_list: ", smile_list)
        # Create a figure with two subplots
        fig, axs = plt.subplots(len(smile_list) + 1, 1, figsize=(8, 6))
        x_axis_max = max(experimental_smile_intensity_dict[smile][0])

        # Iterate over the data and plot each subplot
        for i, (mz_values, intensity_values) in enumerate(smile_list):
            ax = axs[i]  # Select the appropriate subplot

             # Normalize the intensity values
            if len(intensity_values) == 0:
                print("This is an empty intensity_values")
                continue
            max_intensity = np.max(intensity_values)
            normalized_intensity = intensity_values / max_intensity
            #print("normalized_intensity: ", sorted(normalized_intensity))

            # Calculate the width of each bar
#            bar_width = np.min(np.diff(mz_values)) * 800
            bar_width = 0.2

            # Create the bar graph
            ax.bar(mz_values, normalized_intensity, width=bar_width)

            # Set plot labels and title
            ax.set_xlabel('m/z')
            ax.set_ylabel('Intensity')
            ax.set_title(f'm/z Intensity Bar Graph for ' + smile + f'_{i+1}')

            ax.set_xlim([0, x_axis_max + 5])

            # print(" describing_dict[smile]: ",  describing_dict[smile])
            first_elements = [t[0] for t in describing_dict[smile][i]]
            red_indices = first_elements  # Change this to the desired indices

            # this is the elements in the experimenal data
            second_elements = [t[1] for t in describing_dict[smile][i]]

            # Iterate over the red indices and set the color of each corresponding bar to red
            for index in red_indices:
                ax.patches[index].set_color('red')

            if smile not in proteomer_intensity_dict:
               # print("experimental_smile_intensity_dict[smile]", experimental_smile_intensity_dict[smile], second_elements)
                proteomer_intensity_dict[smile] = [second_elements]
            else:
                proteomer_intensity_dict[smile].append([second_elements])
             
        # then plot the experimental result 
        ax = axs[len(smile_list)]
        smile_list = experimental_smile_intensity_dict[smile]
        max_intensity = np.max(smile_list[1])
        normalized_intensity = smile_list[1] / max_intensity

        # Calculate the width of each bar
        #bar_width = np.min(np.diff(smile_list[0])) * 10
        bar_width = 0.2

        # Create the bar graph
        ax.bar(smile_list[0], normalized_intensity, width=bar_width)

        # Set plot labels and title
        ax.set_xlabel('m/z')
        ax.set_ylabel('Intensity')
        ax.set_title(f'experimental result of m/z Intensity Bar Graph for ' + smile)

        unique_second_elements = set()
        [unique_second_elements.add(element[1]) for sublist in describing_dict[smile] for element in sublist]

        unique_second_elements = list(unique_second_elements)

        print("unique_second_elements:", unique_second_elements)
        green_indices = unique_second_elements

        for index in green_indices:
                ax.patches[index].set_color('green')

        ax.set_xlim([0, x_axis_max+5])

        # Adjust spacing between subplots
        fig.tight_layout()       

        # Save the figure
        try:
            plt.savefig(str(smile)+'_intensity.png')
        except:
            print("weird naming problem when saving intensity ")
        plt.clf()

    # plot the barplot
    print("proteomer_intensity_dict", proteomer_intensity_dict)
    


    return experimental_smile_intensity_dict

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_json_smile_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('json_smile_peak_file', type=str, help='json_smile_peak_file is the converted file that has the information that each smiles in the json file with its peak and intensity')
    parser.add_argument('input_mgf_file', type=str, help='input_mgf_file is the mgf file from nextflow')
    parser.add_argument('real_energy_file', type=str, help='input_mgf_file is the mgf file from HPCC')
    parser.add_argument('experiment_for_each_mol', type=str, help='input_mgf_file is the mgf file from HPCC')
    parser.add_argument('tolerance', type=str, help='input_mgf_file is the mgf file from HPCC')

    args = parser.parse_args()

    input_json_smile_file = args.input_json_smile_file
    json_smile_peak_file = args.json_smile_peak_file 
    input_mgf_file = args.input_mgf_file
    real_energy_file = args.real_energy_file
    experiment_for_each_mol = int(args.experiment_for_each_mol)
    tolerance = int(args.tolerance)
    n_peak = 100
    delete_precursor_switch = True

    mgf_dict = get_mgf_dict(input_mgf_file, delete_precursor_switch)
    # print("mgf_dict: ", mgf_dict)

    smile_peak_dict = get_smile_peak_dict(json_smile_peak_file)
    #print("smile_peak_dict: ", smile_peak_dict)
    #return 0
    largest_n_peaks_spect_dict = get_n_largest_peaks(smile_peak_dict, n_peak)

    scan_smile_dict, energy_list, molecule_list = get_scan_smile_dict(input_json_smile_file)
    #print("energy_list: ", energy_list)

    id_energy_dict = get_id_energy_dict(real_energy_file)
    #print("id_energy_dit: ", id_energy_dit)

    smile_id_dict = get_smile_id_dict(json_smile_peak_file)
    #print("smile_id_dict: ", smile_id_dict)

    # mgf_dict:                             scan# <=> mgf's peak_intensity
    # scan_smile_dict: scan# <=> smile        |
    # smile_peak_dict:                      smile <=> json's peak_intensity
    #smile_id_dict:    smile <=> id           |
    # id_energy_dit:                         id   <=> energy
   

    score_list=[]
    energy_score_list = []
    for i in range(0, experiment_for_each_mol):
        energy_score_list.append([])
    missing_counter = 0
    smile_groundTruth_energy_dict = {}

    mgf_keys = []
    for i in mgf_dict.keys():
        # print("i: ", i)
        mgf_keys.append(i.split("-")[0])
    # print("mgf_keys:", mgf_keys)
    # print("mgf_dict.keys()", mgf_dict.keys())

    fig, ax1 = plt.subplots()
    #analysis n_largest_intensity match 
    best_protonation_match_rate = []
    matched_peaks = 0
    molecule_list = []
    describing_dict = {}
    predicted_mgf_dict = {}
    describing_true_peak_dict = {}
    describing_mz_dict = {}

    for i in range(1, len(scan_smile_dict)+1) :
        # if i > 4:
        #     break
        molecule_list.append(i)
        smile = scan_smile_dict[str(i)]
        print("smile: ", smile)
        
        protonation_i = 1
        prediction_lists = []
        combined_experiment_list = []
        
        if str(i) in mgf_keys:        
            counter = 1

            while(True):
                #print("smile_peak_dict: ", smile_peak_dict)
                #return 0
                if smile+"_" + str(counter) not in smile_peak_dict:
                    print("bad", smile+"_" + str(counter))
                    continue

                protonation_index = str(i) + "-" + str(protonation_i)
                print("protonation_index checking",  protonation_index)

                if protonation_index not in mgf_dict.keys():
                    print("protonation_index not in protonation_index",  protonation_index)
                    break

                # manual 
                spec1=mgf_dict[protonation_index]
                #print("spec1: ", spec1)
                spec2=smile_peak_dict[smile+"_" + str(counter)]
                #print(spec2)
                total_experiment_peaks = len(spec2.mz)
                n_peak = total_experiment_peaks
                # top n largest peaks
                #spec2=largest_n_peaks_spect_dict[smile+"_" + str(counter)]
                
                #print("spec2: ", spec2)
                score, matched_peaks=_cosine_fast(spec1,spec2,0.02,False)
                print("score: ", score)


                n_peak_over_1 = len(matched_peaks) / n_peak * 100
                best_protonation_match_rate.append(n_peak_over_1)
                
                print(f"matching peaks out of {n_peak}: {len(matched_peaks)}, maching percenation is: {n_peak_over_1}%", "matched_peaks: ",matched_peaks)
        
                if smile not in describing_dict.keys():
                    describing_dict[smile] = [matched_peaks]
                    predicted_mgf_dict[smile] = [mgf_dict[protonation_index]]
                else:
                    describing_dict[smile].append(matched_peaks)
                    predicted_mgf_dict[smile].append(mgf_dict[protonation_index])

                protonation_i = protonation_i + 1

                intensities = [t[1] for t in matched_peaks]
                if smile not in describing_true_peak_dict.keys():
                    describing_true_peak_dict[smile] = intensities
                else:
                    for j in intensities:
                        if j not in describing_true_peak_dict[smile]:
                            describing_true_peak_dict[smile].append(j)

                with open("table.csv", 'a', newline='') as csvfile:
                    # Create a CSV writer object
                    writer = csv.writer(csvfile)
                    data = [smile, score, len(matched_peaks), matched_peaks],
                    # Write the data to the CSV file
                    writer.writerows(data)

                # print("smile: ", smile)
                # print("check score: ", score)
                # print("len(matched_peaks): ", len(matched_peaks))
                # print("counting smiels: ", smile)
                combined_prediction_list = list(zip(spec1.mz, spec1.intensity))
                combined_experiment_list = list(zip(spec2.mz, spec2.intensity))
                # print(protonation_index, "combined_prediction_list: ", combined_prediction_list)
                # print(protonation_index, "combined_experiment_list: ", combined_experiment_list)
                prediction_lists.append(np.array(combined_prediction_list))


            # ------ calculate weights
            # just checking the first one first, will delete this later:
            
            # #print("prediction_lists: ", prediction_lists)
            # if len(prediction_lists) == 0:
            #     continue
            # all_list = prediction_lists
            # all_list.append(np.array(combined_experiment_list))


            # # padding
            # #lists = [np.pad(lst, ((0, max_length - len(lst)), (0, 0)), mode='constant') for lst in all_list]
            # print("index i = ", i)
            # #print("lists: ", all_list)
            # try:
            #     x_values = np.unique(np.concatenate([lst[:, 0] for lst in all_list ]))
            # except:
            #     print("wrong dimension error")
            #     continue
            # padded_lists = []

            # # Pad and fill the lists with their respective X-axis values
            # for lst in all_list:
            #     padded_lst = np.zeros((len(x_values), 2))
            #     for j, x in enumerate(x_values):
            #         if x in lst[:, 0]:
            #             index = np.where(lst[:, 0] == x)[0][0]
            #             padded_lst[j] = lst[index]
            #     for j, x in enumerate(x_values):
            #         if padded_lst[j, 0] == 0:
            #             padded_lst[j, 0] = x
            #     padded_lists.append(padded_lst)

            # lists = padded_lists

            # initial_weights = np.ones(len(lists) - 1) / (len(lists) - 1)
            # problem = minimize(cosine_score_function, initial_weights, args=(lists,), method='Nelder-Mead')
            # optimized_weights = problem.x

            # # Print the optimized weights and result list
            # print("optimal ratio")
            # for i, weight in enumerate(optimized_weights):
            #     print(f"Ratio for List{i+1}: {weight:.2f}")

            # result_list = np.sum([w * lst for w, lst in zip(optimized_weights, lists[:-1])], axis=0)
            # # print("Result List:")
            # # print(result_list)

            # # Compute and print the cosine score of each list compared to list 5
            # experiment_mz = [x[0] for x in  lists[-1]]
            # experiment_intensity = [x[1] for x in  lists[-1]]
            # experiment_spec = SpectrumTuple("0", experiment_mz, norm_intensity(experiment_intensity))

            # cosine_scores = []
            # for i, lst in enumerate(lists[:-1]):
            #     #get correct format
            #     predict_mz = [x[0] for x in lst]
            #     predict_intensity = [x[1] for x in lst]
            #     predict_spec = SpectrumTuple("0", predict_mz, norm_intensity(predict_intensity))
                                            

            #     cosine_score = cal_cosine_score(predict_spec, experiment_spec, 0.02, True)
            #     cosine_scores.append(cosine_score)

            # print("Cosine Scores:")
            # for i, score in enumerate(cosine_scores):
            #     print(f"List{i+1} vs List5: {score:.4f}")

            # # Compute and print the cosine score between the combined list and list 5
            # combined_list = np.sum([w * np.array(lst) for w, lst in zip(optimized_weights, lists[:-1])], axis=0)

            # predict_mz = [x[0] for x in combined_list]
            # predict_intensity = [x[1] for x in combined_list]
            # predict_spec = SpectrumTuple("0", predict_mz, norm_intensity(predict_intensity))
            # combined_cosine_score = cal_cosine_score(predict_spec, experiment_spec, 0.02, True)
            # print("Combined List vs List5:")
            # print(f"{combined_cosine_score:.4f}")
            

            

    #---- temp:
    describing_intensity_dict = {}
    counter = 1
    ranked_dict = rank_by_intensity(smile_peak_dict, n_peak)
    for smile in describing_true_peak_dict.keys():       
        #smile = "C1=CC=C2C(=C1)C(C(=CC2=O)C)=O"
        describing_mz_dict[smile] = []
        describing_intensity_dict[smile] = []
        print("smile:", smile )
        print(describing_true_peak_dict[smile])
        for i in describing_true_peak_dict[smile]:
            describing_mz_dict[smile].append(smile_peak_dict[smile+"_" + str(1)].mz[i])
        for i in describing_true_peak_dict[smile]:
            describing_intensity_dict[smile].append(smile_peak_dict[smile+"_" + str(1)].intensity[i])
        print("describing_intensity_dict[smile]", describing_intensity_dict[smile])

        # --- non-overlap
        non_overlap_intensity_list = []
        non_overlap_mz_list = []
        overlap_index_set = set(describing_true_peak_dict[smile])
        for index, i in enumerate(smile_peak_dict[smile+"_1"].intensity):
            if index not in overlap_index_set:
                non_overlap_intensity_list.append(i)
                non_overlap_mz_list.append(smile_peak_dict[smile+"_1"].mz[index])
        print("describing_mz_dict[smile]", describing_mz_dict[smile])
        print("len(non_overlap_intensity_list):", len(non_overlap_intensity_list), "non_overlap_intensity_list:", non_overlap_intensity_list)
        print("len(non_overlap_mz_list):", len(non_overlap_mz_list), "non_overlap_mz_list:", non_overlap_mz_list)
        print("all ranked mz", ranked_dict[smile+"_1"].mz)
  
        overlap_intensity = describing_intensity_dict[smile]
        non_overlap_intensity = non_overlap_intensity_list
        bin_range = [0,1]
        number_of_bins = 20

        plt.hist(overlap_intensity, bins = number_of_bins, range = bin_range, edgecolor = 'black', label = "overlap")
        plt.hist(non_overlap_intensity, bins = number_of_bins,range = bin_range,  edgecolor = 'black', alpha = 0.5, label = "non-overlap")

        plt.xlabel("intensity")
        plt.ylabel("frequency")
        plt.legend()

        plt.savefig('hist' + str(counter) + '.png')
        counter += 1
        plt.clf()
            



        # print("describing_dict: ", describing_dict)
        # print("predicted_mgf_dict: ", predicted_mgf_dict)
        

#    return 0
    print("describing_true_peak_dict: ", describing_true_peak_dict)
    for i in describing_true_peak_dict.keys():
        print("len describing_true_peak_dict i: ", i, len(describing_true_peak_dict[i]), len(smile_peak_dict[i+"_1"].mz))
        #print("len smile_peak_dict[smile+ + str(counter)] i: ", i, len(smile_peak_dict[i+"_1"].mz))

    experimental_smile_intensity_dict = plot_discription_function(describing_dict, predicted_mgf_dict, smile_peak_dict)
    print("get here describing_dict")
    plot_barplot(describing_dict, experimental_smile_intensity_dict)
    #plot_discription_rate(ax1, best_protonation_match_rate, molecule_list, "-", "b", "Best_protonation", "upper left")

    print("---------------- Done best protonation ------------")

    #----------- all protonation ------------
    protonation_number = []
    matched_peaks = 0
    total_match_peak = []
    for i in range(1, len(scan_smile_dict)+1) :
        total_match_peak.append([])
        #print(mgf_dict.keys())
        smile = scan_smile_dict[str(i)]
        protonation_i = 1
        count = 0
        if str(i) in mgf_keys:        
            counter = 1
            while(True):
                if smile+"_" + str(counter) not in smile_peak_dict:
                    continue
                protonation_index = str(i) + "-" + str(protonation_i)
                if protonation_index not in mgf_dict.keys():
                    break
                spec1=mgf_dict[protonation_index]
                spec2=largest_n_peaks_spect_dict[smile+"_" + str(counter)]
                print("spec2 2:", spec2)
                score, matched_peaks=_cosine_fast(spec1,spec2,0.02,False)

                n_peak_over_1 = len(matched_peaks) / n_peak * 100
                print("protonation_index: ",protonation_index,f"matching peaks out of {n_peak}: {len(matched_peaks)}, maching percenation is: {n_peak_over_1}%", "matched_peaks: ",matched_peaks)
                protonation_i = protonation_i + 1

                for match in matched_peaks:
                    acutal_peaks = [x for (a,x) in total_match_peak[i-1]]
                    if match[1] not  in acutal_peaks :
                        print("matched not in match_peak: ", match)
                        total_match_peak[i-1].append(match)
                count = count +  1
            protonation_number.append(count)

    for i in range(0, len(total_match_peak)):
        print(i, total_match_peak[i])

    total_match_peak_rate = [len(x)/n_peak*100 for x in total_match_peak]
    print("total_match_peak_rate: ", total_match_peak_rate)
    plot_discription_rate(ax1, total_match_peak_rate, molecule_list, "-", "r", "All_possible_protonation", "upper left")

    print("---------------- Done all protonation ------------")
    ax2 = ax1.twinx()
    plot_discription_rate(ax2, protonation_number, molecule_list,"--", "g", "protonation_times", "upper right")
    plt.savefig('discrption_rate.png')


    


    molecule_energy_score_list = np.array(sorted_list).T
    print(molecule_energy_score_list)
    max_molecule_energy_score_list = []
    max_score_list = []
    useful_molecule_index = []
    print("max_molecule_energy_score_list: ", max_molecule_energy_score_list)
    counter = 0
    for i in molecule_energy_score_list:
        max_score_list.append(max(i))
        max_molecule_energy_score_list.append(energy_list[np.argmax(i)])
        #store the molecule that has the highest score > shreshold
        if max(i) > 0.7:
            useful_molecule_index.append(counter)
            print("useful molecule: ", molecule_list[counter])
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
