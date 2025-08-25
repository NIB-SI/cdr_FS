import pandas as pd
import numpy as np
from scipy.stats import wasserstein_distance
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
standardized_file_rel = config['feature_selection']['standardized_file']
emd_scores_dir_rel = config['feature_selection']['emd_scores_dir']

# Construct full paths
input_file = os.path.join(base_path, standardized_file_rel)
output_dir = os.path.join(base_path, emd_scores_dir_rel)
output_file = os.path.join(output_dir, "EMD_c11_2.5_97.5_well.txt")

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Set trimming percentage
lower_percentile = 2.5
upper_percentile = 97.5

# Read the data
data = pd.read_csv(input_file, sep="\t")

# Check data types
print("Data types:\n", data.dtypes)

# Check unique values in key columns
print("Unique values in Metadata_Day:", data["Metadata_Day"].unique())
print("Unique values in Metadata_Biorep:", data["Metadata_Biorep"].unique())
print("Unique values in Concentration:", data["Concentration"].unique())

# Strip spaces and normalize case in the key columns
data["Metadata_Day"] = data["Metadata_Day"].astype(str).str.strip().str.upper()
data["Metadata_Biorep"] = data["Metadata_Biorep"].str.strip().str.upper()
data["Concentration"] = data["Concentration"].astype(str).str.strip()

# Define the metadata columns and the concentration of interest
metadata_columns = [
    "Concentration", "counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei", 
    "Metadata_Well", "Metadata_Day", "Metadata_Biorep", "Tech_replica", "Day_Well_BR", "cell_ID"
]
concentration = "11"
days = ["D1", "D5", "D7", "D9"]

# Define the combinations for each group, now including BR4
bioreps = ["BR1", "BR2", "BR3", "BR4"]

# Prepare a list to store the results
results_list = []

# Function to trim the extreme values
def trim_extremes(values, lower_percentile=0, upper_percentile=100):
    if lower_percentile == 0 and upper_percentile == 100:
        return values
    if len(values) == 0:
        return values
    lower_bound = np.percentile(values, lower_percentile)
    upper_bound = np.percentile(values, upper_percentile)
    return values[(values >= lower_bound) & (values <= upper_bound)]

# Iterate over each day
for day in days:
    print(f"\nProcessing day: {day}")
    for i in range(4):  # Now iterate over 4 bioreps
        biorep = bioreps[i]
        print(f"Filter for popu{i+1}: Metadata_Biorep == {biorep}, Metadata_Day == {day}, Concentration == {concentration}")
        filtered_data = data[
            (data["Metadata_Biorep"] == biorep) &
            (data["Metadata_Day"] == day) &
            (data["Concentration"] == concentration)
        ]
        print(f"Filtered data (popu{i+1}):\n", filtered_data.head())

    # Create the population dataframes for the current day and concentration
    populations = {
        f"popu{i+1}": data[
            (data["Metadata_Biorep"] == bioreps[i]) &
            (data["Metadata_Day"] == day) &
            (data["Concentration"] == concentration)
        ] for i in range(4)  # Now create 4 populations
    }

    # Trim the data at the well level
    for key in populations:
        pop = populations[key]
        trimmed_pop = pd.DataFrame()
        for well in pop["Metadata_Well"].unique():
            well_data = pop[pop["Metadata_Well"] == well]
            trimmed_well_data = well_data.copy()
            for feature in well_data.columns:
                if feature not in metadata_columns:
                    trimmed_well_data[feature] = trim_extremes(well_data[feature].dropna(), lower_percentile, upper_percentile)
            trimmed_pop = pd.concat([trimmed_pop, trimmed_well_data])
        populations[key] = trimmed_pop

    for key, pop in populations.items():
        print(f"{key} size after trimming: {pop.shape[0]} rows")
        # Print sample of trimmed data
        if pop.shape[0] > 0:
            print(pop.head())

    # Iterate over each pair of populations
    for i in range(4):  # Now iterate over 4 populations
        for j in range(i+1, 4):  # And compare with the remaining populations
            popu1 = populations[f"popu{i+1}"]
            popu2 = populations[f"popu{j+1}"]
            population1 = f"{bioreps[i]}_{day}_{concentration}"
            population2 = f"{bioreps[j]}_{day}_{concentration}"

            print(f"Comparing: {population1} vs {population2}")

            # Iterate over each feature column (ignoring metadata columns)
            for feature in data.columns:
                if feature not in metadata_columns:
                    values1 = popu1[feature].dropna()
                    values2 = popu2[feature].dropna()

                    # Check if both distributions are completely missing
                    if values1.empty and values2.empty:
                        print(f"Skipping feature {feature} - no data in both populations")
                        continue

                    # Remove NaN values that may appear after trimming
                    values1 = values1.dropna()
                    values2 = values2.dropna()

                    # Check if both distributions are completely missing after trimming
                    if values1.empty or values2.empty:
                        print(f"Skipping feature {feature} after trimming - no data in one or both populations")
                        continue

                    # Calculate EMD ignoring missing values
                    emd_score = wasserstein_distance(values1, values2)
                    results_list.append({
                        "Feature": feature,
                        "Population1": population1,
                        "Population2": population2,
                        "Concentration": concentration,
                        "Metadata_Day": day,
                        "Count1": len(values1),
                        "Count2": len(values2),
                        "EMD_score": emd_score
                    })
                    print(f"Feature: {feature}, EMD: {emd_score}, Count1: {len(values1)}, Count2: {len(values2)}")

# Convert the results list to a DataFrame
results = pd.DataFrame(results_list)

# Save the results to a tab-separated file
results.to_csv(output_file, sep="\t", index=False)

print("EMD controls analysis completed and saved to", output_file)