import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import configparser
from scipy.optimize import curve_fit

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

# Get paths from config
base_path = config['feature_selection']['base_path']
emd_scores_dir_rel = config['feature_selection']['emd_scores_dir']
emd_scores_drc_dir_rel = config['feature_selection']['emd_scores_drc_dir']

# Construct full paths
file_path = os.path.join(base_path, emd_scores_dir_rel, 'EMD_conc_2.5_97.5_well.txt')
output_dir = os.path.join(base_path, emd_scores_drc_dir_rel)

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Step 2: Read the data
data = pd.read_csv(file_path, sep='\t')

# Step 3: Define the population pairs and order
population_pairs = [
    (11, 10),
    (11, 9),
    (11, 8),
    (11, 7),
    (11, 6),
    (11, 5),
    (11, 4),
    (11, 3),
]
pair_labels = [f'{p1}v{p2}' for p1, p2 in population_pairs]

# Step 4: Extract unique features and days
features = data['Feature'].unique()
days = data['Metadata_Day'].unique()

# Step 5: Define the models for dose-response fitting
def brain_cousens_bc4(x, b, c, d, e):
    return c + (d - c) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))

def brain_cousens_bc5(x, b, c, d, e, f):
    return c + (d - c + f * x) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))

def four_param_log_logistic(x, b, c, d, e):
    return c + (d - c) / (1 + np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))

def four_param_weibull(x, b, c, d, e):
    return c + (d - c) * np.exp(-np.exp(b * (np.log(x + 1e-10) - np.log(e + 1e-10))))

def linear(x, m, b):
    return m * x + b

def constant(x, c):
    return np.full_like(x, c, dtype=np.float64)

# Step 6: Fit model and calculate AIC and BIC
def fit_dose_response(x, y, initial_guess, model='BC4'):
    # Remove NaN values for fitting
    mask = ~np.isnan(y)
    x = x[mask]
    y = y[mask]

    if len(x) > 0:
        try:
            # Fit the model
            if model == 'BC4':
                params, cov = curve_fit(brain_cousens_bc4, x, y, p0=initial_guess, maxfev=50000)
                residuals = y - brain_cousens_bc4(x, *params)
            elif model == 'BC5':
                params, cov = curve_fit(brain_cousens_bc5, x, y, p0=initial_guess, maxfev=50000)
                residuals = y - brain_cousens_bc5(x, *params)
            elif model == 'LL4':
                params, cov = curve_fit(four_param_log_logistic, x, y, p0=initial_guess, maxfev=50000)
                residuals = y - four_param_log_logistic(x, *params)
            elif model == 'WB1.4':
                params, cov = curve_fit(four_param_weibull, x, y, p0=initial_guess, maxfev=50000)
                residuals = y - four_param_weibull(x, *params)
            elif model == 'Lin':
                params, cov = curve_fit(linear, x, y, p0=initial_guess, maxfev=50000)
                residuals = y - linear(x, *params)
                slope = params[0]  # Slope (k value)
            elif model == 'Con':
                params, cov = curve_fit(constant, x, y, p0=[np.mean(y)], maxfev=50000)
                residuals = y - constant(x, *params)
            else:
                raise ValueError("Unknown model type specified.")

            # Calculate AIC and BIC
            ss_res = np.sum(residuals**2)
            k = len(params)
            n = len(y)
            aic = n * np.log(ss_res / n) + 2 * k
            bic = n * np.log(ss_res / n) + k * np.log(n)
            
            if model == 'Lin':
                return params, aic, bic, slope  # Return slope for linear model
            else:
                return params, aic, bic, None
        except (RuntimeError, OverflowError) as e:
            print(f"Warning: Optimal parameters not found for x={x}, y={y}. Error: {str(e)}")
            return None, None, None, None
    else:
        return None, None, None, None

# Step 7: Plotting for one feature and one day
def plot_feature_day(feature, day, data, results):
    print(f"Processing feature: {feature}, day: {day}")  # Debug statement
    subset = data[(data['Feature'] == feature) & (data['Metadata_Day'] == day)]
    x = np.arange(len(population_pairs))
    y = []

    for p1, p2 in population_pairs:
        score = subset[((subset['Population1'] == p1) & (subset['Population2'] == p2)) |
                       ((subset['Population1'] == p2) & (subset['Population2'] == p1))]['EMD_score']
        if not score.empty:
            y.append(score.values[0])
        else:
            y.append(np.nan)

    y = np.array(y, dtype=np.float64)
    plt.plot(x, y, 'o')

    models = ['BC4', 'BC5', 'LL4', 'WB1.4', 'Lin', 'Con']
    colors = ['#d73027', '#f46d43', '#fdae61', '#abd9e9', '#74add1', '#4575b4']
    initial_guesses = [
        [1, min(y), max(y), np.median(x)],
        [1, min(y), max(y), np.median(x), 0.1],
        [1, min(y), max(y), np.median(x)],
        [1, min(y), max(y), np.median(x)],
        [1, np.mean(y)],
        [np.mean(y)]
    ]

    legend_entries = []

    # Check for NaNs before fitting the models
    if not np.any(np.isnan(y)):
        for model, color, initial_guess in zip(models, colors, initial_guesses):
            params, aic, bic, slope = fit_dose_response(x, y, initial_guess, model)
            if params is not None:
                if model == 'Lin':
                    results.append((feature, day, model, aic, bic, slope))
                    plt.plot(x, linear(x, *params), '-', color=color, alpha=0.9)
                    legend_entries.append((model, aic, bic, color, slope))
                else:
                    results.append((feature, day, model, aic, bic, None))
                    if model == 'BC4':
                        plt.plot(x, brain_cousens_bc4(x, *params), '-', color=color, alpha=0.9)
                    elif model == 'BC5':
                        plt.plot(x, brain_cousens_bc5(x, *params), '-', color=color, alpha=0.9)
                    elif model == 'LL4':
                        plt.plot(x, four_param_log_logistic(x, *params), '-', color=color, alpha=0.9)
                    elif model == 'WB1.4':
                        plt.plot(x, four_param_weibull(x, *params), '-', color=color, alpha=0.9)
                    elif model == 'Con':
                        plt.plot(x, constant(x, *params), '-', color=color, alpha=0.9)
                    legend_entries.append((model, aic, bic, color, None))

    plt.title(f'{feature}, Day: {day}', fontsize=8, wrap=True)
    plt.ylabel('EMD Score', fontsize=8)
    plt.xticks(ticks=x, labels=pair_labels, rotation=45, fontsize=6)
    
    # Sort legend entries by AIC + BIC
    legend_entries.sort(key=lambda x: x[1] + x[2])
    handles = [plt.Line2D([0], [0], color=color, lw=2, alpha=0.9) for _, _, _, color, _ in legend_entries]
    labels = [
        f'{model} AIC/BIC: {aic:.2f}/{bic:.2f}' if model != 'Lin' else f'{model} AIC/BIC: {aic:.2f}/{bic:.2f}, Slope: {slope:.2f}'
        for model, aic, bic, color, slope in legend_entries
    ]
    plt.legend(handles, labels, fontsize=8)

# Step 8: Generate plots for all features and days and save results
def generate_all_plots(features, days, data, output_dir):
    results = []
    for day in days:
        num_features = len(features)
        num_plots = 25  # 5x5 grid
        num_images = (num_features + num_plots - 1) // num_plots
        
        for img_idx in range(num_images):
            plt.figure(figsize=(20, 20))
            for i in range(num_plots):
                feature_idx = img_idx * num_plots + i
                if feature_idx < num_features:
                    feature = features[feature_idx]
                    plt.subplot(5, 5, i + 1)
                    plot_feature_day(feature, day, data, results)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'Day_{day}_part_{img_idx + 1}.png'))
            plt.close()
    
    # Save results to a tab-separated .txt file
    results_df = pd.DataFrame(results, columns=['Feature', 'Day', 'Model', 'AIC', 'BIC', 'Slope'])
    results_df.to_csv(os.path.join(output_dir, 'model_fit_results.txt'), sep='\t', index=False)

# Call the function to generate all plots
generate_all_plots(features, days, data, output_dir)

print("DRC model fitting and plotting completed. Results saved to", output_dir)