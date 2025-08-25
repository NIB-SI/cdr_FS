import pandas as pd
import numpy as np
from scipy.stats import wasserstein_distance
import os
import configparser
from itertools import combinations

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
trimmed_standardized_file_rel = config['feature_selection']['trimmed_standardized_file']
emd_scores_dir_rel = config['feature_selection']['emd_scores_dir']

# Construct full paths
input_file = os.path.join(base_path, trimmed_standardized_file_rel)
output_dir = os.path.join(base_path, emd_scores_dir_rel)
output_file = os.path.join(output_dir, "EMD_conc_2.5_97.5_well.txt")

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Set trimming percentage
lower_percentile = 0
upper_percentile = 100

# Read the data
data = pd.read_csv(input_file, sep="\t")

# Strip spaces and normalize case in the key columns
data["Metadata_Day"] = data["Metadata_Day"].astype(str).str.strip().str.upper()
data["Metadata_Biorep"] = data["Metadata_Biorep"].str.strip().str.upper()
data["Concentration"] = data["Concentration"].astype(str).str.strip()

# Define the metadata columns
metadata_columns = [
    "Concentration", "counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei", 
    "Metadata_Well", "Metadata_Day", "Metadata_Biorep", "Tech_replica", "Day_Well_BR", "cell_ID"
]

# List of concentrations to compare
concentrations = ["2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]
days = ["D1", "D5", "D7", "D9"]

# Prepare a list to store the results
results_list = []

# Function to trim the extreme values
def trim_extremes(values, lower_percentile=0, upper_percentile=100):
    if lower_percentile == 0 and upper_percentile == 100:
        return values
    if len(values) == 0:
        return values
    lower_bound = np.nanpercentile(values, lower_percentile)
    upper_bound = np.nanpercentile(values, upper_percentile)
    return values[(values >= lower_bound) & (values <= upper_bound)]

# Function to compute EMD and store results
def compute_emd(day, conc1, conc2, data):
    features = data.columns.difference(metadata_columns)
    for feature in features:
        values1 = data[(data["Metadata_Day"] == day) & (data["Concentration"] == conc1)][feature]
        values2 = data[(data["Metadata_Day"] == day) & (data["Concentration"] == conc2)][feature]

        # Check if both distributions are completely missing
        if values1.empty and values2.empty:
            print(f"Skipping feature {feature} - no data in both populations")
            continue

        # Remove NaN values
        values1 = values1.dropna()
        values2 = values2.dropna()

        # Trim extreme values
        values1 = trim_extremes(values1, lower_percentile, upper_percentile)
        values2 = trim_extremes(values2, lower_percentile, upper_percentile)

        # Check if both distributions are completely missing after trimming
        if values1.empty or values2.empty:
            print(f"Skipping feature {feature} after trimming - no data in one or both populations")
            continue

        emd_score = wasserstein_distance(values1, values2)
        results_list.append({
            "Feature": feature,
            "Population1": conc1,
            "Population2": conc2,
            "Concentration": f"{conc1}_vs_{conc2}",
            "Metadata_Day": day,
            "Count1": len(values1),
            "Count2": len(values2),
            "EMD_score": emd_score
        })
        print(f"Feature: {feature}, EMD: {emd_score}, Count1: {len(values1)}, Count2: {len(values2)}")

# Generate all pairwise comparisons of concentrations for each day
for day in days:
    for conc1, conc2 in combinations(concentrations, 2):
        compute_emd(day, conc1, conc2, data)

# Create a DataFrame from the results and save it to a file
results_df = pd.DataFrame(results_list)
results_df.to_csv(output_file, sep="\t", index=False)

print("EMD calculations completed and saved to", output_file)