import scanpy as sc
import numpy as np
import os
import warnings

warnings.filterwarnings('ignore')

# Exact EBV Gene List from GTF
EXACT_EBV_GENES = [
    'LMP-2B', 'BNRF1', 'EBER-1_pol_III_transcript', 'EBER-2_pol_III_transcript', 'BCRF1', 
    'BWRF1', 'EBNA-2', 'BHLF1', 'BHRF1', 'BFLF2', 'BFLF1', 'BFRF1A', 'BFRF1', 'BFRF2', 
    'BFRF3', 'BPLF1', 'BOLF1', 'BORF1', 'BORF2', 'BaRF1', 'BMRF1', 'BMRF2', 'BSLF2-BMLF1', 
    'BSLF1', 'BSRF1', 'BLLF3', 'BLRF1', 'BLRF2', 'BLLF1', 'BLLF2', 'EBNA-3A', 'EBNA-3B-EBNA-3C', 
    'BZLF2', 'BZLF1', 'BRLF1', 'BRRF1', 'BRRF2', 'EBNA-1', 'BKRF2', 'BKRF3', 'BKRF4', 
    'BBLF4', 'BBRF1', 'BBRF2', 'BBLF2-BBLF3', 'BBRF3', 'BBLF1', 'BGLF5', 'BGLF4', 'BGLF3.5', 
    'BGLF3', 'BGRF1-BDRF1', 'BGLF2', 'BGLF1', 'BDLF4', 'BDLF3.5', 'BDLF3', 'BDLF2', 'BDLF1', 
    'BcLF1', 'BcRF1', 'BTRF1', 'BXLF2', 'BXLF1', 'BXRF1', 'BVRF1', 'BVLF1', 'BVRF2', 'BdRF1', 
    'BILF2', 'RPMS1', 'LF3', 'LF2', 'LF1', 'BILF1', 'BALF5', 'BALF4', 'A73', 'BALF3', 'BARF0', 
    'BALF2', 'BALF1', 'BARF1', 'BNLF2b', 'BNLF2a', 'LMP-1'
]

def clean_and_recluster(filepath, output_suffix="_post_qc.h5ad"):
    print(f"\n{'='*50}\nProcessing: {os.path.basename(filepath)}\n{'='*50}")
    
    adata = sc.read_h5ad(filepath)
    initial_cells = adata.n_obs
    print(f"Initial cells: {initial_cells:,}")

    # 1. Calculate QC metrics if not present
    if 'n_genes_by_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    n_genes_col = 'n_genes_by_counts' if 'n_genes_by_counts' in adata.obs.columns else 'n_genes'

    # 2. Hard Gene Filters (Leena's Request)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=6000)
    print(f"Cells after gene thresholds (200-6000): {adata.n_obs:,}")

    # 3. Scrublet Doublet Removal
    print("Running Scrublet to identify doublets...")
    sc.external.pp.scrublet(adata, expected_doublet_rate=0.05)
    adata = adata[adata.obs['predicted_doublet'] == False].copy()
    print(f"Cells after Doublet Removal: {adata.n_obs:,}")

    # 4. Strict Viral UMI >= 2 Filter
    print("Calculating strict GTF viral load...")
    adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
    adata.var_names_make_unique()
    
    ebv_genes = [g for g in adata.var_names if g in EXACT_EBV_GENES or g.startswith('EBV-')]
    
    if adata.raw is not None:
        raw_var_names = adata.raw.var_names.str.replace(r'(?i)_type1$', '', regex=True)
        ebv_raw_indices = [i for i, g in enumerate(raw_var_names) if g in ebv_genes]
        raw_viral_sum = np.ravel(adata.raw.X[:, ebv_raw_indices].sum(axis=1))
    else:
        raw_viral_sum = np.ravel(adata[:, ebv_genes].X.sum(axis=1))

    # Apply UMI >= 2 logic: If it has 1 UMI, set count to 0 (treat as background/ambient)
    adata.obs['Viral_Counts'] = np.where(raw_viral_sum >= 2, raw_viral_sum, 0)
    
    if 'Infection' in adata.obs.columns:
        adata.obs['Infection'] = np.where(adata.obs['Viral_Counts'] > 0, 'EBV', 'Mock')

    print(f"True EBV+ Cells (UMI >= 2): {np.sum(adata.obs['Viral_Counts'] > 0):,}")

    # 5. Enhanced Clustering & Separation (Leena's Request)
    print("Recalculating PCA and UMAP with high-separation parameters...")
    
    # We must operate on the highly variable genes for PCA
    if 'highly_variable' not in adata.var.columns:
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        
    sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True, n_comps=50)
    sc.pp.neighbors(adata, n_pcs=50, n_neighbors=15)
    
    # min_dist=0.3 and spread=1.5 will forcefully separate B and T cells visually
    sc.tl.umap(adata, min_dist=0.3, spread=1.5)

    # 6. Save Output
    out_path = filepath.replace(".h5ad", output_suffix)
    adata.write_h5ad(out_path)
    print(f"✅ Cleaned file saved to: {out_path}")
    return out_path

# ==========================================
# Run the Pipeline
# ==========================================
FILE_WT = "/Users/schung/work/app_sc_MS_WT_NC007605/optimized_data/data_200_annotated.h5ad"
FILE_ENR = "/Users/schung/work/comp_MS_WT_vs_EN_NC007606/optimized_data/data_EBV_enriched_100.h5ad"

new_wt = clean_and_recluster(FILE_WT)
new_enr = clean_and_recluster(FILE_ENR)

print("\n🎉 Post-QC Cleanup Complete!")
print(f"Update your integration script to use:\n1. {new_wt}\n2. {new_enr}")