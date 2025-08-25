import os
import configparser

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
correlation_dir_rel = config['feature_selection']['correlation_dir']
correlation_list_include_dir_rel = config['feature_selection']['correlation_list_include_dir']

# Construct full paths
input_dir = os.path.join(base_path, correlation_dir_rel)
output_dir = os.path.join(base_path, correlation_list_include_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List all input files in the input directory
input_files = [f for f in os.listdir(input_dir) if f.endswith('.txt')]

# Process each input file
for input_file in input_files:
    # Define the input file path
    input_file_path = os.path.join(input_dir, input_file)
    
    # Read the input file and extract the first feature from each cluster
    clean_features = []
    with open(input_file_path, 'r') as f:
        for line in f:
            if line.strip():
                cluster, features = line.split(':')
                features_list = eval(features.strip())
                if features_list:
                    clean_features.append(features_list[0])
    
    # Define the output file path
    output_file_path = os.path.join(output_dir, input_file.replace("high_corr_clusters", "clean_features"))
    
    # Save the clean features to the output file
    with open(output_file_path, 'w') as f:
        for feature in clean_features:
            f.write(f'{feature}\n')
    
    print(f"Processed and saved clean features for {input_file} to {output_file_path}")

print(f"All clean features saved to {output_dir}")