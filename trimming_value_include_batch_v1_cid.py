import pandas as pd
import numpy as np
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
standardized_file_rel = config['feature_selection']['standardized_file']
trimmed_dir_rel = config['feature_selection']['trimmed_dir']
emd_scores_drc_selected_dir_rel = config['feature_selection']['emd_scores_drc_selected_dir']

# Construct full paths
input_file = os.path.join(base_path, standardized_file_rel)
output_dir = os.path.join(base_path, trimmed_dir_rel)
feature_lists_dir = os.path.join(base_path, emd_scores_drc_selected_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List all feature list files in the feature lists directory
feature_lists_files = [f for f in os.listdir(feature_lists_dir) if f.endswith('.txt')]

# Set trimming percentage
lower_percentile = 2.5
upper_percentile = 97.5

# Read the input data
data = pd.read_csv(input_file, sep="\t")

# Strip spaces and normalize case in the key columns
data["Metadata_Day"] = data["Metadata_Day"].astype(str).str.strip().str.upper()
data["Metadata_Biorep"] = data["Metadata_Biorep"].str.strip().str.upper()
data["Concentration"] = data["Concentration"].astype(str).str.strip()
data["Metadata_Well"] = data["Metadata_Well"].astype(str).str.strip().str.upper()
data["Day_Well_BR"] = data["Day_Well_BR"].astype(str).str.strip().str.upper()
data["cell_ID"] = data["cell_ID"].astype(str).str.strip().str.upper()

# Define the metadata columns
metadata_columns = [
    "Concentration", "counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei", 
    "Metadata_Well", "Metadata_Day", "Metadata_Biorep", "Tech_replica", "Day_Well_BR", "cell_ID"
]

# Function to trim the extreme values
def trim_extremes(values, lower_percentile=0, upper_percentile=100):
    if len(values) == 0:
        return values
    lower_bound = np.nanpercentile(values, lower_percentile)
    upper_bound = np.nanpercentile(values, upper_percentile)
    return values[(values >= lower_bound) & (values <= upper_bound)]

# Process each feature list
for feature_list_file in feature_lists_files:
    # Load the inclusion list
    include_features_path = os.path.join(feature_lists_dir, feature_list_file)
    include_features = pd.read_csv(include_features_path, header=None)[0].tolist()

    # Include only the features from the main data
    included_features = metadata_columns + [col for col in data.columns if col in include_features]
    data_filtered = data[included_features]

    # Prepare a list to store the results
    trimmed_data_list = []

    # Trim the data at the well level
    grouped = data_filtered.groupby(["Metadata_Day", "Metadata_Biorep", "Metadata_Well"])
    for (day, biorep, well), group in grouped:
        trimmed_group = group.copy()
        for feature in group.columns.difference(metadata_columns):
            trimmed_values = trim_extremes(group[feature].dropna(), lower_percentile, upper_percentile)
            trimmed_group[feature] = trimmed_values
        trimmed_data_list.append(trimmed_group)

    # Concatenate all trimmed data into a single DataFrame
    trimmed_data_well = pd.concat(trimmed_data_list)

    # Define the output file path
    output_file_well = os.path.join(output_dir, feature_list_file.replace("features", "trimmed_features"))

    # Save the trimmed data to a file
    trimmed_data_well.to_csv(output_file_well, sep="\t", index=False)

    print(f"Trimming by Well level completed and saved to {output_file_well}")

print(f"All trimmed datasets saved to {output_dir}")