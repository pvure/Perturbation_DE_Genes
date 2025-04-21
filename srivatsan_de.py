import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re 


processed_adata_path = "/Users/pranayvure/Srivatsan_DE/srivatsan20_processed.h5ad" 
import anndata
import pandas as pd
import numpy as np
# No need for 're' if we only check for a literal comma

# Define the path to the processed AnnData file
processed_adata_path = "srivatsan20_processed.h5ad"

try:
    print(f"Attempting to load data from: {processed_adata_path}")
    adata = anndata.read_h5ad(processed_adata_path)
    print(f"Successfully loaded data: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- Part 1: Identify and Show Examples Containing Commas ---
    print("\n--- Identifying Perturbations Containing Commas ---")

    if 'perturbation' in adata.obs.columns:
        separator_to_check = ','

        # Use regex=False for simple string checking
        comma_mask = adata.obs['perturbation'].astype(str).str.contains(separator_to_check, regex=False)
        
        num_comma_cells = comma_mask.sum()

        if num_comma_cells > 0:
            print(f"Found {num_comma_cells} cells where 'perturbation' contains a comma ('{separator_to_check}').")
            
            print("\nA few unique examples:")
            comma_obs = adata.obs[comma_mask]
            unique_examples = comma_obs['perturbation'].unique()
            print(unique_examples[:min(len(unique_examples), 10)]) # Show up to 10 unique examples
        else:
             print(f"Found 0 cells where 'perturbation' contains a comma ('{separator_to_check}').")

    else:
         print("\nSkipped investigation: 'perturbation' column not found in adata.obs.")

    # --- Part 2: Code to Remove Cells with Commas ---
    print("\n--- Code Snippet for Removing Cells with Commas ---")
    if 'perturbation' in adata.obs.columns and num_comma_cells > 0:

        # Create a boolean mask: True for cells TO REMOVE (containing comma)
        # We recalculate it here just to be clear in the final snippet
        mask_to_remove_comma = adata.obs['perturbation'].astype(str).str.contains(',', regex=False)

        # Calculate how many cells will be kept vs removed
        cells_to_remove_count = mask_to_remove_comma.sum()
        cells_to_keep_count = (~mask_to_remove_comma).sum()

        print(f"\nOriginal number of cells: {adata.n_obs}")
        print(f"Number of cells to be removed (contain comma): {cells_to_remove_count}")
        print(f"Number of cells to be kept (do NOT contain comma): {cells_to_keep_count}")

        print("\nCode to perform the filtering (removes cells with commas):")
        print("```python")
        print("# --- Filtering Code (Comma Removal) ---")
        print("# Ensure 'adata' is your loaded AnnData object")
        print("")
        print("# Create boolean mask (True for cells to REMOVE - those containing a comma)")
        print("mask_to_remove_comma = adata.obs['perturbation'].astype(str).str.contains(',', regex=False)")
        print("")
        print("# Keep cells that do NOT contain a comma (using ~)")
        print("adata_filtered = adata[~mask_to_remove_comma, :].copy()")
        print("")
        print("# Verify the new shape")
        print("print(f'Original shape: {adata.shape}')")
        print("print(f'Filtered shape: {adata_filtered.shape}')")
        print("```")

    elif 'perturbation' not in adata.obs.columns:
        print("\nCannot provide filtering code: 'perturbation' column not found.")
    elif num_comma_cells == 0:
        print("\nNo cells found containing commas. No filtering needed for this condition.")


except FileNotFoundError:
    print(f"Error: The file '{processed_adata_path}' was not found.")
    print("Please ensure the processed file exists in the correct location.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    # import traceback
    # traceback.print_exc()