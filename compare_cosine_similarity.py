import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pyteomics import mgf
import collections
from typing import List, Tuple
import itertools
import argparse

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)

def norm_intensity(intensity):
    return np.copy(intensity)/np.linalg.norm(intensity)

def _cosine_fast(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    allow_shift: bool,
) -> Tuple[float, List[Tuple[int, int]]]:
    precursor_charge = max(spec.precursor_charge, 1)
    precursor_mass_diff = (
        spec.precursor_mz - spec_other.precursor_mz
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

def main():
    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('input_structures_file', type=str, help='input_structures_file')
    parser.add_argument('output_file', type=str, help='output_file')

    args = parser.parse_args()

    input_structures_file = args.input_structures_file
    output_file = args.output_file

    delete_precursor_switch = False
    spec_dic = {}
    count = 0
    scan_max = 0
    for spectrum in mgf.read(input_structures_file):
        print(spectrum)
        params = spectrum.get('params')
        mz = spectrum.get('m/z array')
        intensity = spectrum.get('intensity array')
        spec_dic[params['scan']] = {}

        scan_max = max(scan_max, int(params['scan']))
        if not delete_precursor_switch:
            spec_dic[params['scan']] = SpectrumTuple(params['pepmass'][0],int(params['charge'][0]),mz,norm_intensity(intensity))
        else: 
            precursor = params['pepmass'][0]
            out_mz_intensity = []

            mz_out = []
            intensity_out = []
            for i in range(0, len(mz)):
                if mz[i] < float(precursor) - 5:
                    mz_out.append(mz[i])
                    intensity_out.append(intensity[i])

            print("np.array(mz_out): ", mz_out)
            spec_dic[params['scan']] = SpectrumTuple(params['pepmass'][0],int(params['charge'][0]), mz_out, norm_intensity(intensity_out))
        #print(spec_dic[params['scan']])
        count = count + 1

    score_list = []
    score_per_mol = []
    print("len(spec_dic): " , len(spec_dic))
    for i in range(0,int(scan_max/4)):
        score_per_mol.append([])
        score = 0
        for comb in itertools.combinations([x+1 for x in range(4)],2):
            print(str(comb[0]))
            try:
                spec1=spec_dic[str(comb[0] + i*4)]
                spec2=spec_dic[str(comb[1] + i*4)]
                score, matched_peaks=_cosine_fast(spec1,spec2,0.1,False)
                score_list.append(score)
                score_per_mol[i].append(score)
            except:
                continue

    print(score_per_mol)

    fig,ax=plt.subplots()
    ax.boxplot(score_per_mol)

    ax.set_title("results")
    ax.set_xlabel("molecules")
    ax.set_ylabel("consine similarity * 100%")

    plt.savefig(output_file + '.png')


if __name__ == "__main__":
    main()