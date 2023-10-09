# this is the readme for how to use these scripts. (Yunshu wrote these scripts for himself to use.) 

# Convert MSn mgf to Json file:
```
python ~/QCXMS/scripts/convert_mgf_to_json.py 20230811_mce_library_pos_all_lib_MSn.mgf 20230811_mce_library_pos_all_lib_MSn.json

```

# Merge multiple proteomer result to a single molecule:
```
python ~/QCXMS/scripts/merge_multi_proteomers.py merged_predicted_spectra.mgf merged_multi_proteomers.txt
```