import anndata
import pandas as pd
import numpy as np
import os # Import os module to check if file exists

# --- Configuration ---
# Input file path (the curated but not yet filtered file)
input_adata_path = "srivatsan20_processed.h5ad"
# Output file path (where the filtered data will be saved)
output_adata_path = "srivatsan20_filtered.h5ad"
# Separator indicating multiple perturbations to remove
separator_to_remove = ','

# --- Main Logic ---
try:
    # --- Load Data ---
    print(f"Attempting to load data from: {input_adata_path}")
    if not os.path.exists(input_adata_path):
         raise FileNotFoundError(f"Input file not found: {input_adata_path}")
    adata = anndata.read_h5ad(input_adata_path)
    print(f"Successfully loaded data: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- Identify Cells to Remove ---
    if 'perturbation' not in adata.obs.columns:
        raise KeyError("Column 'perturbation' not found in adata.obs. Cannot perform filtering.")

    print(f"\nIdentifying cells where 'perturbation' contains '{separator_to_remove}'...")
    # Create boolean mask (True for cells to REMOVE - those containing the separator)
    mask_to_remove = adata.obs['perturbation'].astype(str).str.contains(separator_to_remove, regex=False)
    num_to_remove = mask_to_remove.sum()

    if num_to_remove == 0:
        print(f"Found 0 cells containing '{separator_to_remove}'. No filtering needed based on this criterion.")
        print(f"Saving the original data to '{output_adata_path}' as no filtering was applied.")
        # If no filtering is needed, you might still want to save it under the new name
        # or skip saving. Here, we save it anyway for consistency.
        adata_filtered = adata.copy() # Create a copy to avoid modifying original if needed elsewhere
    else:
        print(f"Found {num_to_remove} cells containing '{separator_to_remove}'. These will be removed.")

        # --- Filter Data ---
        print("\nFiltering data...")
        # Keep cells that do NOT contain the separator (using ~)
        # Use .copy() to ensure adata_filtered is a new object in memory
        adata_filtered = adata[~mask_to_remove, :].copy()
        print(f"Original shape: {adata.shape}")
        print(f"Filtered shape: {adata_filtered.shape}")

    # --- Save Filtered Data ---
    print(f"\nSaving filtered data to: {output_adata_path}")
    # Add compression to potentially reduce file size
    adata_filtered.write(output_adata_path, compression='gzip')
    print("Filtered data saved successfully.")

except FileNotFoundError as fnf_error:
    print(f"Error: {fnf_error}")
    print("Please ensure the input file exists at the specified path.")
except KeyError as key_error:
     print(f"Error: {key_error}")
     print("Please check the structure of your input AnnData object.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    # import traceback
    # traceback.print_exc()