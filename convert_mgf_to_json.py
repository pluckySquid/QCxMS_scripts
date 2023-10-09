# useage
import json
import argparse

def parse_mgf(mgf_text):
    """Parse an MGF formatted text and return a list of entries."""
    entries = mgf_text.strip().split("BEGIN IONS")[1:]
    parsed_entries = []
    for entry in entries:
        entry_data = {}
        peaks = []
        lines = entry.strip().split("\n")
        for line in lines:
            if "=" in line:
                key, value = line.split("=", 1)
                entry_data[key] = value
            elif "END IONS" not in line and len(line.split()) == 2:
                mz, intensity = line.split()
                peaks.append([float(mz), float(intensity)])
        entry_data["PEAKS"] = peaks
        parsed_entries.append(entry_data)
    return parsed_entries

def convert_to_json(entry, spectrum_id):
    """Convert a parsed MGF entry into the given JSON format."""
    return {
        "spectrum_id": spectrum_id,
        "source_file": "N/A",
        "task": "N/A",
        "scan": "N/A",
        "ms_level": entry.get("MSLEVEL", "N/A"),
        "library_membership": "N/A",
        "spectrum_status": "N/A",
        "peaks_json": json.dumps(entry.get("PEAKS", [])),
        "Compound_Name": entry.get("NAME", "N/A"),
        "PI": entry.get("PI", "N/A"),
        "Data_Collector": entry.get("DATACOLLECTOR", "N/A"),
        "Adduct": entry.get("ADDUCT", "N/A"),
        "Precursor_MZ": entry.get("PEPMASS", "N/A"),
        "Charge": entry.get("CHARGE", "N/A"),
        "Smiles": entry.get("SMILES", "N/A"),
        "INCHI": entry.get("INCHI", "N/A"),
        "INCHI_AUX": entry.get("INCHIAUX", "N/A"),
        "SpectrumID": "MSn" + str(spectrum_id),
        "Ion_Mode": entry.get("IONMODE", "N/A"),
        "description": entry.get("DESCRIPTION", "N/A"),
        "exact_mass": entry.get("EXACTMASS", "N/A"),
        "formula": entry.get("FORMULA", "N/A"),
        "feature_id": entry.get("FEATURE_ID", "N/A"),
        "retention_time_in_seconds": entry.get("RTINSECONDS", "N/A"),
        "spec_type": entry.get("SPECTYPE", "N/A"),
        "collision_energy": entry.get("Collision energy", "N/A"),
        "fragmentation_method": entry.get("FRAGMENTATION_METHOD", "N/A"),
        "isolation_window": entry.get("ISOLATION_WINDOW", "N/A"),
        "acquisition": entry.get("Acquisition", "N/A"),
        "instrument_type": entry.get("INSTRUMENT_TYPE", "N/A"),
        "source_instrument": entry.get("SOURCE_INSTRUMENT", "N/A"),
        "ion_source": entry.get("ION_SOURCE", "N/A"),
        "dataset_id": entry.get("DATASET_ID", "N/A"),
        "usi": entry.get("USI", "N/A"),
        "scans": entry.get("SCANS", "N/A"),
        "precursor_purity": entry.get("PRECURSOR_PURITY", "N/A"),
        "quality_chimeric": entry.get("QUALITY_CHIMERIC", "N/A"),
        "quality_explained_intensity": entry.get("QUALITY_EXPLAINED_INTENSITY", "N/A"),
        "quality_explained_signals": entry.get("QUALITY_EXPLAINED_SIGNALS", "N/A")
    }

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Merge all the MGF Files and create a summary CSV')
    parser.add_argument('mgf_file', type=str, help='mgf_file1')
    parser.add_argument('json_file', type=str, help='mgf_file1')


    args = parser.parse_args()

    mgf_file = args.mgf_file
    json_file = args.json_file

    # Load the MGF file
    with open(mgf_file, "r") as f:
        mgf_text = f.read()

    parsed_entries = parse_mgf(mgf_text)
    json_entries = []
    spectrum_id = 1
    for entry in parsed_entries:
        json_entries.append(convert_to_json(entry, spectrum_id))
        spectrum_id += 1

    # Write to the JSON file
    with open(json_file, "w") as f:
        json.dump(json_entries, f, indent=4)

    print(f"Conversion complete! JSON data has been saved to {json_file}.")
