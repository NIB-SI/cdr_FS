import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import configparser
from scipy.cluster.hierarchy import dendrogram, linkage, set_link_color_palette
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
trimmed_dir_rel = config['feature_selection']['trimmed_dir']
correlation_dir_rel = config['feature_selection']['correlation_dir']

# Construct full paths
input_dir = os.path.join(base_path, trimmed_dir_rel)
output_dir = os.path.join(base_path, correlation_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# List all input files in the input directory
input_files = [f for f in os.listdir(input_dir) if f.endswith('.txt')]

# Function to trim extreme values
def trim_extremes(values, lower_percentile=0, upper_percentile=100):
    if len(values) == 0:
        return values
    lower_bound = np.nanpercentile(values, lower_percentile)
    upper_bound = np.nanpercentile(values, upper_percentile)
    return values[(values >= lower_bound) & (values <= upper_bound)]

# Process each input file
for input_file in input_files:
    # Define the input file path
    input_file_path = os.path.join(input_dir, input_file)
    
    # Load data
    data = pd.read_csv(input_file_path, sep='\t')
    
    # Normalize case and strip spaces in key columns
    data["Metadata_Day"] = data["Metadata_Day"].astype(str).str.strip().str.upper()
    data["Metadata_Biorep"] = data["Metadata_Biorep"].str.strip().str.upper()
    data["Concentration"] = data["Concentration"].astype(str).str.strip()
    data["Metadata_Well"] = data["Metadata_Well"].astype(str).str.strip().str.upper()
    data["Day_Well_BR"] = data["Day_Well_BR"].astype(str).str.strip().str.upper()
    data["cell_ID"] = data["cell_ID"].astype(str).str.strip().str.upper()
    
    # Create the "Well_level" column
    data["Well_level"] = data["Metadata_Day"] + "_" + data["Metadata_Biorep"] + "_" + data["Metadata_Well"]
    
    # Define the metadata columns, including the new "Well_level" column
    metadata_columns = [
        "Concentration", "counts_Cells", "counts_Cytoplasm", "counts_FilteredNuclei", 
        "Metadata_Well", "Metadata_Day", "Metadata_Biorep", "Tech_replica", "Well_level", "Day_Well_BR", "cell_ID"
    ]
    
    # Fill NaN values with the mean of each column for numeric data
    data = data.apply(lambda x: x.fillna(x.mean()) if x.dtype != 'object' else x)
    
    # Select only non-metadata columns for median calculation
    numeric_columns = data.columns.difference(metadata_columns)
    
    # Group by "Well_level" and calculate the median values for numeric columns only
    grouped_data = data.groupby("Well_level")[numeric_columns].median().reset_index()
    
    # Calculate Pearson correlation matrix, excluding "Well_level"
    feature_columns = grouped_data.columns.difference(["Well_level"])
    features = grouped_data[feature_columns]
    correlation_matrix = features.corr(method='pearson')
    
    # Handle missing and infinite values in the distance matrix
    dist_matrix = 1 - np.abs(correlation_matrix)
    dist_matrix = np.nan_to_num(dist_matrix, nan=1.0, posinf=1.0, neginf=1.0)  # Replace NaN and inf with 1.0
    
    # Ensure the distance matrix contains only valid values for hierarchical clustering
    dist_matrix[dist_matrix < 0] = 0  # Set any negative values to 0
    condensed_dist_matrix = squareform(dist_matrix, checks=False)
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(condensed_dist_matrix, method='median')
    
    # Use AgglomerativeClustering to form clusters based on the correlation threshold
    agg_cluster = AgglomerativeClustering(n_clusters=None, metric='precomputed', linkage='average', distance_threshold=0.1)
    clusters = agg_cluster.fit_predict(dist_matrix)
    
    # Combine features into clusters
    feature_clusters = pd.DataFrame({'Feature': feature_columns, 'Cluster': clusters})
    high_corr_clusters = feature_clusters.groupby('Cluster')['Feature'].apply(list)
    
    # Save the clusters to a file
    output_txt_file = os.path.join(output_dir, f'high_corr_clusters_{input_file.replace(".txt", "")}.txt')
    with open(output_txt_file, 'w') as f:
        for cluster, features in high_corr_clusters.items():
            f.write(f'Cluster {cluster}: {features}\n')
    
    # Create a color palette for the clusters with more than two members
    significant_clusters = high_corr_clusters[high_corr_clusters.apply(len) > 2]
    unique_significant_clusters = significant_clusters.index
    colors = sns.color_palette('husl', len(unique_significant_clusters))
    hex_colors = sns.color_palette('husl', len(unique_significant_clusters)).as_hex()
    cluster_colors = {cluster: color for cluster, color in zip(unique_significant_clusters, hex_colors)}
    
    # Set the color palette for the dendrogram
    set_link_color_palette([cluster_colors[cluster] for cluster in unique_significant_clusters])
    
    # Default color for single or two-member clusters
    default_color = 'gray'
    
    # Create a dictionary to map each feature to its cluster color, using default for small clusters
    feature_to_color = {}
    for cluster, features in significant_clusters.items():
        color = cluster_colors[cluster]
        for feature in features:
            feature_to_color[feature] = color
    
    # Use default color for features not in significant clusters
    for feature in feature_clusters['Feature']:
        if feature not in feature_to_color:
            feature_to_color[feature] = default_color
    
    # Extract the list of feature names
    features_list = feature_columns.tolist()
    
    # Create a dendrogram with increased size and resolution, and color the labels and lines
    plt.figure(figsize=(25, 15), dpi=300)
    dendro = dendrogram(linkage_matrix, labels=features_list, leaf_rotation=90, leaf_font_size=10, above_threshold_color='gray', color_threshold=0.1)
    
    # Color the labels based on the clusters
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(feature_to_color.get(lbl.get_text(), default_color))
    
    plt.title('Hierarchical Clustering Dendrogram Well Level')
    plt.xlabel('Feature')
    plt.ylabel('Distance')
    plt.axhline(y=0.1, color='r', linestyle='--')
    plt.xticks(rotation=90, ha='right', fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'dendrogram_median_well_{input_file.replace(".txt", "")}.png'))
    plt.close()
    
    print(f"Results for {input_file} saved to {output_dir}")

print(f"All correlation analysis results saved to {output_dir}")