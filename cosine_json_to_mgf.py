# in ~/QCXMS/data/ 
# python ~/QCXMS/scripts/cosine_json_to_mgf.py mols_GNPS-SELLECKCHEM-FDA-PART1_first_350.txt GNPS-SELLECKCHEM-FDA-PART1_smile_peaks merged_predicted_spectra.mgf
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
import numpy as np


SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["mz", "intensity"]
)

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
                peak_match_scores.append(
                    peak_intensity * spec_other.intensity[other_peak_i]
                )
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


def get_mgf_dict(input_mgf_file):
    spec_dic = {}
    count = 0
    for spectrum in mgf.read(input_mgf_file):
        #print(spectrum)
        params = spectrum.get('params')
        mz = spectrum.get('m/z array')
        intensity = spectrum.get('intensity array')
        spec_dic[params['scan']] = {}
        spec_dic[params['scan']] = SpectrumTuple(mz,norm_intensity(intensity))
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
                data = line.split(": ")[1].split(": ")[0]
                
                pattern = r'\d+\.\d+'
                matches = re.findall(pattern, data)

                # Convert the matches to float and group them into pairs
                coords = [[float(matches[i]), float(matches[i+1])] for i in range(0, len(matches), 2)]

                mz = [x[0] for x in coords]
                #print(mz)
                intensity = [x[1] for x in coords]
                #print(intensity)

                smile_peak_dic[smile] = SpectrumTuple(mz,norm_intensity(intensity))
    #print(smile_peak_dic)
    return smile_peak_dic

def get_scan_smile_dict(input_json_smile_file): 
    scan_smile_dict = {}
    with open(input_json_smile_file) as input_json_smile_file:
        scanNumber = 0
        for line in input_json_smile_file:
            if line != "smiles,energy,atom\n":
                scan_smile_dict[str(scanNumber)] = line.split(",")[0]
            scanNumber = scanNumber + 1

    return scan_smile_dict

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_json_smile_file', type=str, help='input_json_smile_file is the txt file that has the information about which all the inpt smiles')
    parser.add_argument('json_smile_peak_file', type=str, help='json_smile_peak_file is the converted file that has the information that each smiles in the json file with its peak and intensity')
    parser.add_argument('input_mgf_file', type=str, help='input_mgf_file is the mgf file from HPCC')
    

    args = parser.parse_args()

    input_json_smile_file = args.input_json_smile_file
    json_smile_peak_file = args.json_smile_peak_file 
    input_mgf_file = args.input_mgf_file

    mgf_dict = get_mgf_dict(input_mgf_file)
    #print(mgf_dict)

    smile_peak_dict = get_smile_peak_dict(json_smile_peak_file)
    #print(smile_peak_dic)

    scan_smile_dict = get_scan_smile_dict(input_json_smile_file)
    #print(scan_smile_dict)

    score_list=[]
    energy_score_list = [[],[],[],[],[],[],[]]
    energy_counter  = 1
    for key, values in mgf_dict.items():
        tmp_score = []
        
        smile = scan_smile_dict[str(key)]
        #print(smile)
        counter = 1
        while True:
            if smile+"_" + str(counter) not in smile_peak_dict:
                #print(smile+"_" + str(counter))
                #print(smile_peak_dict[smile+"_" + str(counter)])
                break
            
            spec1=values
            spec2=smile_peak_dict[smile+"_" + str(counter)]
            score, matched_peaks=_cosine_fast(spec1,spec2,0.1,False)
            tmp_score.append(score)
            counter = counter + 1
            #print(tmp_score)

            energy_score_list[energy_counter % 7].append(score)

        print(str(key) + ": " + str(tmp_score))
        score_list.append(tmp_score)
        energy_counter += 1

    print("energy_counter: ", energy_counter)
    ave_list = []
    sorted_list = []
    for i in range(0, 7):
        #sum = sum(energy_score_list[i])
        ave_list.append(sum(energy_score_list[i]) / len(energy_score_list[i]))
        sorted_list.append(sorted(energy_score_list[i]))
       
        #print(sorted_list)
    # note the ave_list's index zer0 is enery = 70
    print("length are ")
    print(len(x) for x in energy_score_list)
    print(sorted_list)

    print("len(energy_score_list[i])", len(energy_score_list[i]))

    mol_energy_list = np.array(energy_score_list).T

    colors = ['red', 'blue', 'green', 'purple', 'orange', 'brown', 'pink']
    for i in range(len(sorted_list)):
        if i != 0:
            plt.plot(energy_score_list[i], color=colors[i], label=f"energy = {i * 10}")
        else:
            plt.plot(energy_score_list[i], color=colors[i], label="energy = 70")

    # set plot title and labels
    plt.title('Seven Lists in One Plot')
    plt.xlabel('# number experiments')
    plt.ylabel('cosine')
    plt.legend()
    plt.savefig('sorted_list.png')  # Save the plot as a PNG file

# add legend
plt.legend()



if __name__ == "__main__":
    main()