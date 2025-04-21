import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define the path to the FILTERED AnnData file
filtered_adata_path = "/Users/pranayvure/Srivatsan_DE/srivatsan20_filtered.h5ad" # Make sure this file exists

try:
    print(f"Attempting to load filtered data from: {filtered_adata_path}")
    adata_filtered = anndata.read_h5ad(filtered_adata_path)
    print(f"Successfully loaded filtered data: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} genes")

    # --- Analyze Dose Frequencies per Perturbation ---
    print("\n--- Analyzing Dose Frequencies ---")

    if 'perturbation' in adata_filtered.obs.columns and 'dose_value' in adata_filtered.obs.columns:
        # Work with the .obs dataframe
        obs_df = adata_filtered.obs.copy()

        # Filter out control cells ('vehicle') if they still exist
        # *** Corrected line below using .str.lower() ***
        obs_df = obs_df[obs_df['perturbation'].astype(str).str.lower() != 'vehicle']

        if obs_df.empty:
            print("No chemically perturbed cells found in the filtered data.")
        else:
            # Ensure 'dose_value' is numeric, coercing errors to NaN
            obs_df['dose_value_numeric'] = pd.to_numeric(obs_df['dose_value'], errors='coerce')

            # Drop rows where dose could not be converted (optional, but good practice)
            num_nan_doses = obs_df['dose_value_numeric'].isna().sum()
            if num_nan_doses > 0:
                 print(f"Warning: {num_nan_doses} rows had non-numeric dose values and were excluded.")
                 obs_df = obs_df.dropna(subset=['dose_value_numeric'])

            if obs_df.empty:
                 print("No valid numeric dose values found after cleaning.")
            else:
                # --- 1. Frequency Table ---
                print("\n--- Generating Dose Frequency Table ---")
                # Group by perturbation and numeric dose value, count cells
                dose_counts = obs_df.groupby(['perturbation', 'dose_value_numeric']).size()

                # Unstack to create a table: perturbations (rows) vs doses (columns)
                # fill_value=0 ensures doses not used for a specific perturbation show 0 count
                dose_freq_table = dose_counts.unstack(fill_value=0)

                # Sort columns (doses) for better readability
                dose_freq_table = dose_freq_table.reindex(sorted(dose_freq_table.columns), axis=1)

                print("Frequency Table (Perturbations vs Dose Values, showing first 10 perturbations):")
                print(dose_freq_table.head(10).to_string()) # Use to_string to potentially show more columns

                # Optional: Identify the most common dose for each perturbation
                most_common_doses = dose_freq_table.idxmax(axis=1)
                print("\nMost common dose per perturbation (first 10):")
                print(most_common_doses.head(10))


                # --- 2. Heatmap Visualization ---
                print("\n--- Generating Dose Frequency Heatmap ---")

                # Check if the table is too large for effective heatmap visualization
                n_perturbations = dose_freq_table.shape[0]
                n_doses = dose_freq_table.shape[1]
                print(f"Heatmap will have {n_perturbations} perturbations (rows) and {n_doses} doses (columns).")

                # Decide whether to plot all or subset (e.g., if > 50 perturbations)
                if n_perturbations > 50:
                    print("Warning: Large number of perturbations. Heatmap may be crowded.")
                    # Optional: Subset the table, e.g., to top N perturbations by total cell count
                    # top_n = 50
                    # top_pert_indices = dose_freq_table.sum(axis=1).nlargest(top_n).index
                    # plot_table = dose_freq_table.loc[top_pert_indices]
                    # print(f"Plotting top {top_n} perturbations by cell count.")
                    plot_table = dose_freq_table # Keep all for now, user can uncomment/modify subsetting
                else:
                    plot_table = dose_freq_table

                if plot_table.empty:
                     print("Cannot generate heatmap: No data to plot.")
                else:
                    plt.figure(figsize=(max(8, n_doses * 0.6), max(10, n_perturbations * 0.2))) # Adjust size dynamically

                    # Use log(count + 1) for color scale if counts vary widely, helps visualization
                    log_counts = np.log1p(plot_table)

                    sns.heatmap(log_counts, cmap="viridis", annot=False, # Annot=True might be too crowded
                                linewidths=.5, linecolor='lightgray', cbar_kws={'label': 'log(Cell Count + 1)'})

                    plt.title("Frequency of Dose Values per Perturbation")
                    plt.xlabel("Dose Value")
                    plt.ylabel("Perturbation")
                    # Make y-axis labels smaller if too many perturbations
                    if n_perturbations > 50:
                         plt.tick_params(axis='y', labelsize=8)
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout() # Adjust layout
                    plt.show()

except FileNotFoundError:
    print(f"Error: The file '{filtered_adata_path}' was not found.")
    print("Please ensure the filtered file exists.")
except KeyError as key_error:
     print(f"Error: Required column not found - {key_error}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    # import traceback
    # traceback.print_exc()