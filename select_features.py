import pandas as pd
import numpy as np
import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
emd_scores_drc_dir_rel = config['feature_selection']['emd_scores_drc_dir']
emd_scores_drc_selected_dir_rel = config['feature_selection']['emd_scores_drc_selected_dir']

# Construct full paths
file_path = os.path.join(base_path, emd_scores_drc_dir_rel, 'model_fit_results.txt')
output_path = os.path.join(base_path, emd_scores_drc_selected_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Load the data
data = pd.read_csv(file_path, sep="\t")

# Define day combinations
day_combinations = {
    #"D1": ["D1"],
    #"D5": ["D5"],
    #"D7": ["D7"],
    #"D9": ["D9"],
    #"D1_D5_D7": ["D1", "D5", "D7"],
    #"D5_D7": ["D5", "D7"],
    "all_days": ["D1", "D5", "D7", "D9"],
}

# Function to check if a feature should be excluded
def should_exclude_feature(feature_data):
    for day, day_data in feature_data.groupby('Day'):
        day_data['AIC+BIC'] = day_data['AIC'] + day_data['BIC']
        if (day_data.loc[day_data['Model'] == 'Con', 'AIC+BIC'].min() == day_data['AIC+BIC'].min()).any():
            return True
    return False

# Initialize a dictionary to store the filtered data for each day combination
filtered_data = {}

# Filter data for each day combination
for combo_name, days in day_combinations.items():
    combo_data = data[data['Day'].isin(days)]
    
    # Filter out features with non-positive Slope values in the Lin model for the current combination
    valid_features = combo_data[(combo_data['Model'] == 'Lin') & (pd.to_numeric(combo_data['Slope'], errors='coerce') > 0)]['Feature'].unique()
    combo_data_filtered = combo_data[combo_data['Feature'].isin(valid_features)]
    
    # Check exclusion criteria for each feature
    final_valid_features = []
    for feature, feature_data in combo_data_filtered.groupby('Feature'):
        if not should_exclude_feature(feature_data):
            final_valid_features.append(feature)
    
    filtered_data[combo_name] = final_valid_features

# Save the lists for each day combination
for combo_name, features in filtered_data.items():
    output_file_path = os.path.join(output_path, f"{combo_name}_features.txt")
    with open(output_file_path, "w") as file:
        for feature in features:
            file.write(f"{feature}\n")
    print(f"Saved {len(features)} features for {combo_name} to {output_file_path}")

print(f"All feature selection lists saved to {output_path}")