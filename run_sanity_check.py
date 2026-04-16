import scanpy as sc

# 1. Load the new combined file
print("Loading new combined file...")
adata_combined = sc.read_h5ad("/Users/schung/work/comp_MS_WT_vs_EN_NC007606/combined_integrated_analysis_MULTI_compartment.h5ad")

# 2. Count the WT Cycling T cells in the new file
wt_mask = adata_combined.obs['dataset'] == 'WT'
new_wt_counts = adata_combined[wt_mask].obs['majority_voting'].value_counts()
print("\n--- NEW COMBINED FILE (WT ONLY) ---")
print(f"Total WT Cells: {wt_mask.sum():,}")
print(f"Cycling T Cells: {new_wt_counts.get('Cycling T cells', 0):,}")

# 3. Load your original, untouched WT file
print("\nLoading original WT file...")
adata_original = sc.read_h5ad("/Users/schung/work/app_sc_MS_WT_NC007605/optimized_data/data_200_annotated.h5ad", backed='r')

# 4. Count the Cycling T cells in the original file
original_counts = adata_original.obs['majority_voting'].value_counts()
print("\n--- ORIGINAL WT FILE ---")
print(f"Total WT Cells: {adata_original.n_obs:,}")
print(f"Cycling T Cells: {original_counts.get('Cycling T cells', 0):,}")