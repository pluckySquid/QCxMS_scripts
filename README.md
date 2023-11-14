# this is the readme for how to use these scripts. (Yunshu wrote these scripts for himself to use.) 

# Analysis
## To get the box plots:
```
python ~/QCXMS/QCxMS_scripts/analysis_sets/analysis_sets.py mols_ecom_test.txt ~/QCXMS/workDir/neccessary_files/FDA-Part2/delete_precursor_17 merged_predicted_spectra.mgf ~/QCXMS/workDir/neccessary_files/ALL_GNPS_cleaned.csv ~/QCXMS/workDir/neccessary_files/FDA-Part2/GNPS-SELLECKCHEM-FDA-PART2.json 8
```

## Make bar plots
```
python ~/QCXMS/QCxMS_scripts/NIH-Nature_round2_positive/analysis_auto_protonation_heatmap_backup_green_on_exp.py mols_multi_proteomer.txt ~/QCXMS/workDir/neccessary_files/FDA-Part2/delete_precursor_17 merged_predicted_spectra.mgf ~/QCXMS/workDir/neccessary_files/ALL_GNPS_cleaned.csv 1 3
```


## Convert MSn mgf to Json file:
```
python ~/QCXMS/scripts/convert_mgf_to_json.py 20230811_mce_library_pos_all_lib_MSn.mgf 20230811_mce_library_pos_all_lib_MSn.json

```

# Tools

## Merge multiple proteomer result to a single molecule:
```
python ~/QCXMS/scripts/merge_multi_proteomers.py merged_predicted_spectra.mgf merged_multi_proteomers.txt
```

## Get a random mgf from a mgf:
```
python ~/QCXMS/scripts/random_bin/generate_random_peaks.py merged_predicted_spectra.mgf
```

## Generate possible peaks (Reza)
```
python ~/QCXMS/QCxMS_scripts/generate_spectrum/main.py mols_ecom_test.txt predicted.mgf
```

## Convert a crest result to a multi-proteomer result
```
python ~/QCXMS/scripts/convert_crestmgf_to_multiproteomermgf.py merged_predicted_spectra.mgf 
```
## choose top n intensities from a mgf and write to a mgf
```
python ~/QCXMS/scripts/choose_high_intensity_from_mgf.py input.mgf output.mgf
```

## combine crest and multi-proteomers:
```
python ~/QCXMS/scripts/combine_cret_multi_proteomers.py merged_multi_proteomers.txt ../n_time_crest/merged_predicted_spectra.mgf combined_crest_multi.txt
```

## Get a delete_17 file from json file.

```
python ~/QCXMS/scripts/GNPS_json_to_mgf.py ALL_GNPS.json delete_precursor_17
```

# Manually Analysis

## transform csv result, make each line to 8 line and combien them
```
python ~/QCXMS/scripts/manually_analysis/transform_csv.py path_to_input_csv path_to_output_csv
```



# normalize the data:

## find the lowest number of peaks amoung same molecule k, then output each molecule with k peaks. 
This script will first process the data. It will remove mz<50 Da. and mz > precursor - 17. \
Then it will output the processed data to processed_molecules.txt. \
Then it will find the smallest peaks k in each molecule. Let's say, for the first molecule, there are 8 experiments, with peaks: 135
139
141
124
118
153
99
114
This k will be 99.\
Then this will proccess all the experiments for this molecule, that only ouptut the top k intensities to the output file.

```
python ~/QCXMS/QCxMS_scripts/normalization/find_smallest_peaks_for_same_exp.py merged_predicted_spectra.mgf 8 low_peak_spec.mgf
```


## Use the specific peaks to normalize the data (Use this to choose the highest peaks for QCXMS)
```
python ~/QCXMS/QCxMS_scripts/normalization/normalize_all.py merged_predicted_spectra.mgf 8 low_peak_spec.mgf
```


## Use the specific peaks to normalize the data then get the quirtile:
```
python ~/QCXMS/QCxMS_scripts/normalization/quirtile.py merged_predicted_spectra.mgf 8 low_peak_spec.mgf
```

## Use the specific peaks to normalize the data (Use this to choose the highest peaks, normalize_Rezas )
```
python ~/QCXMS/QCxMS_scripts/normalization/normalize_Rezas_prediction.py merged_predicted_spectra.mgf 5 low_peak_spec.mgf
```


# reproducibility
## check 10 molecules samilarity:
```
python ~/QCXMS/QCxMS_scripts/reproducibility/analysis_reproducibility.py \
mols_reproducibility_60_traj.txt \
~/QCXMS/workDir/neccessary_files/FDA-Part2/delete_precursor_17 \
merged_predicted_spectra.mgf \
~/QCXMS/workDir/neccessary_files/ALL_GNPS_cleaned.csv \
~/QCXMS/workDir/neccessary_files/FDA-Part2/GNPS-SELLECKCHEM-FDA-PART2.json \
10
```


## plot the line plot to show trend
```
python ~/QCXMS/QCxMS_scripts/reproducibility/plot_reproducibility_all.py reproducibility.txt reproducibility_all.png
```