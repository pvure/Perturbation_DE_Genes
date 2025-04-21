import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re 


processed_adata_path = "/Users/pranayvure/Srivatsan_DE/srivatsan20_processed.h5ad" 

try:
    print(f"Attempting to load data from: {processed_adata_path}")
    adata = anndata.read_h5ad(processed_adata_path)
    print(f"Successfully loaded data: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # --- Part 1: Show examples for each separator ---
    print("\n--- Investigating Separators in 'perturbation' column ---")
    
    if 'perturbation' in adata.obs.columns:
        # Define the separators we checked previously
        separators = [',', ';', '+', '&', '/']
        
        found_any_examples = False
        for sep in separators:
            # Escape the separator for regex and find cells containing it
            pattern = re.escape(sep)
            matching_obs = adata.obs[adata.obs['perturbation'].astype(str).str.contains(pattern, regex=True)]
            
            if not matching_obs.empty:
                found_any_examples = True
                print(f"\nExamples containing separator '{sep}' (found in {len(matching_obs)} cells):")
                # Get unique perturbation names containing this separator, show up to 5 examples
                unique_examples = matching_obs['perturbation'].unique()
                print(unique_examples[:min(len(unique_examples), 5)]) # Show first 5 unique examples
            #else:
            #    print(f"\nNo perturbation names found containing separator '{sep}'") # Optional: report separators not found

        if not found_any_examples:
             print("\nNo examples found for any specified separator.")

    else:
         print("\nSkipped investigation: 'perturbation' column not found in adata.obs.")

    # --- Part 2: Code to Remove Cells with Separators ---
    print("\n--- Code Snippet for Removing Cells with Separators ---")
    if 'perturbation' in adata.obs.columns and found_any_examples:
        # Define separators and create the combined regex pattern again
        separators_to_remove = [',', ';', '+', '&', '/'] # MODIFY THIS LIST based on investigation
        print(f"The code below will remove cells whose 'perturbation' contains ANY of: {separators_to_remove}")
        print(">>> Please examine the examples above and MODIFY `separators_to_remove` if needed <<<")

        pattern_to_remove = '|'.join(re.escape(sep) for sep in separators_to_remove)

        # Create a boolean mask: True for cells TO REMOVE
        mask_to_remove = adata.obs['perturbation'].astype(str).str.contains(pattern_to_remove, regex=True)

        # Calculate how many cells will be kept vs removed
        cells_to_remove_count = mask_to_remove.sum()
        cells_to_keep_count = (~mask_to_remove).sum()

        print(f"\nOriginal number of cells: {adata.n_obs}")
        print(f"Number of cells to be removed (contain separators): {cells_to_remove_count}")
        print(f"Number of cells to be kept (do NOT contain separators): {cells_to_keep_count}")

        print("\nCode to perform filtering:")
        print("```python")
        print("# --- Filtering Code ---")
        print("# Ensure 'adata' is your loaded AnnData object")
        print("# Ensure 'pattern_to_remove' is defined as above, based on separators you want to filter out")
        print(f"pattern_to_remove = '{pattern_to_remove}' # Based on separators: {separators_to_remove}")
        print("")
        print("# Create boolean mask (True for cells to REMOVE)")
        print("mask_to_remove = adata.obs['perturbation'].astype(str).str.contains(pattern_to_remove, regex=True)")
        print("")
        print("# Keep cells that do NOT match the removal mask (using ~)")
        print("adata_filtered = adata[~mask_to_remove, :].copy()")
        print("")
        print("# Verify the new shape")
        print("print(f'Original shape: {adata.shape}')")
        print("print(f'Filtered shape: {adata_filtered.shape}')")
        print("```")

    elif 'perturbation' not in adata.obs.columns:
        print("\nCannot provide filtering code: 'perturbation' column not found.")
    elif not found_any_examples:
        print("\nNo cells flagged for removal based on current separators. No filtering code needed.")


except FileNotFoundError:
    print(f"Error: The file '{processed_adata_path}' was not found.")
    print("Please ensure the processed file exists in the correct location.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")
    # import traceback
    # traceback.print_exc()