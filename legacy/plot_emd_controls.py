import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import configparser
from matplotlib import gridspec

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
emd_scores_dir_rel = config['feature_selection']['emd_scores_dir']

# Construct full paths
input_file = os.path.join(base_path, emd_scores_dir_rel, 'EMD_c11_2.5_97.5_well.txt')
output_dir = os.path.join(base_path, emd_scores_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

def create_emd_plot(data, feature_order, feature_positions, output_file, figsize=(28.8, 15), is_small=False):
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3], hspace=0.05)

    ax1 = plt.subplot(gs[0])  # Top subplot (log scale)
    ax2 = plt.subplot(gs[1])  # Bottom subplot (linear scale)

    dot_size = 0.2 if is_small else 6
    alpha = 0.6 if is_small else 0.5
    fontsize_factor = 0.8 if is_small else 1.5

    for feature in feature_order:
        feature_data = data[data['Feature'] == feature]
        x = [feature_positions[feature]] * len(feature_data)
        y = feature_data['EMD_score']
        
        mask_top = y > 7
        if np.any(mask_top):
            ax1.scatter(np.array(x)[mask_top], y[mask_top], c='black', s=dot_size, alpha=alpha)
        
        mask_bottom = y <= 7
        ax2.scatter(np.array(x)[mask_bottom], y[mask_bottom], c='black', s=dot_size, alpha=alpha)

    ax1.set_yscale('log')
    ax1.set_ylim(7, data['EMD_score'].max() * 1.1)
    ax1.set_xlim(-1, len(feature_order))
    ax1.set_xticks([])
    
    ax2.set_ylim(0, 7)
    ax2.set_xlim(-1, len(feature_order))

    # Remove only the top and right spines
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if is_small:
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['bottom'].set_linewidth(0.5)

    if not is_small:
        ax2.set_xticks(range(len(feature_order)))
        ax2.set_xticklabels(feature_order, rotation=90, fontsize=4.2*fontsize_factor)
        ax2.tick_params(axis='x', which='major', pad=1)
        ax2.tick_params(axis='y', which='major', labelsize=14*fontsize_factor)
        ax1.tick_params(axis='y', which='major', labelsize=14*fontsize_factor)
    else:
        ax2.set_xticks([])
        ax2.set_xticklabels([])

    # Adjust subplot positions
    if is_small:
        plt.subplots_adjust(left=0.08, top=0.95, bottom=0.1)
    else:
        plt.subplots_adjust(left=0.06, top=0.95, bottom=0.15, right=0.99)  # Adjusted left and right for full-size plot

    # Add y-axis label
    if is_small:
        fig.text(0.02, 0.5, 'EMD score (linear/log scale)', fontsize=14*fontsize_factor, ha='left', va='center', rotation='vertical')
    else:
        fig.text(0.01, 0.5, 'EMD score (linear/log scale)', fontsize=14*fontsize_factor, ha='left', va='center', rotation='vertical')

    # Add 'Features' label
    if is_small:
        fig.text(0.5, 0.02, 'Features', fontsize=12*fontsize_factor, ha='center', va='center')
    else:
        ax2.set_xlabel('Features', fontsize=14*fontsize_factor, labelpad=10)

    fig.suptitle('EMD scores for EMD_c11_2.5_97.5', fontsize=16*fontsize_factor, y=0.98)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Load and process data
data = pd.read_csv(input_file, sep='\t')
feature_order = data.groupby('Feature')['EMD_score'].median().sort_values().index
feature_positions = {feature: i for i, feature in enumerate(feature_order)}

# Create full-size plot
full_size_output = os.path.join(output_dir, 'EMD_c11_2.5_97.5_well_scatter_plot_full_size.png')
create_emd_plot(data, feature_order, feature_positions, full_size_output)

# Create smaller plot for academic article
small_size_output = os.path.join(output_dir, 'EMD_c11_2.5_97.5_well_scatter_plot_small_size.png')
create_emd_plot(data, feature_order, feature_positions, small_size_output, figsize=(8.64, 4.5), is_small=True)

print(f'Full-size EMD scatter plot saved to: {full_size_output}')
print(f'Small-size EMD scatter plot saved to: {small_size_output}')