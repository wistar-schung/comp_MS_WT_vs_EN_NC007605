import scanpy as sc
import numpy as np
import scipy.sparse as sp

print("Loading combined dataset...")
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
adata = sc.read_h5ad(FILE_PATH)

# Clean var names to ensure we catch all viral transcripts
adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
adata.var_names_make_unique()

# Isolate WT only
adata_wt = adata[adata.obs['dataset'] == 'WT'].copy()

# Exact EBV Gene List extracted from EBV_NC_007605_ncRNA_ParseBio.gtf
KNOWN_EBV_GENES = [
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

ebv_genes = [g for g in adata_wt.var_names if g in KNOWN_EBV_GENES or g.startswith('EBV-')]

# --- CRITICAL: Use RAW unnormalized data for read counts if available ---
if adata_wt.raw is not None:
    print("Using raw, unnormalized counts matrix for accurate read percentages...")
    matrix_to_use = adata_wt.raw.X
    # Need to get indices of ebv_genes in the raw var_names
    ebv_raw_indices = [i for i, g in enumerate(adata_wt.raw.var_names.str.replace(r'(?i)_type1$', '', regex=True)) if g in ebv_genes]
    ebv_matrix = matrix_to_use[:, ebv_raw_indices]
else:
    print("WARNING: adata.raw not found. If adata.X is log-normalized, read percentages will be inaccurate.")
    matrix_to_use = adata_wt.X
    ebv_matrix = adata_wt[:, ebv_genes].X

# 1. Total Cells
total_cells = adata_wt.n_obs

# 2. Viral Reads per Cell
viral_reads_per_cell = np.ravel(ebv_matrix.sum(axis=1))

# 3. Infected Cells (Filter: > 0 viral UMIs)
# You can change this to > 1 if you want to be stricter against ambient RNA
infected_cells = np.sum(viral_reads_per_cell >= 2)
pct_ebv_cells = (infected_cells / total_cells) * 100

# 4. Read Calculations
total_viral_reads = np.sum(viral_reads_per_cell)
total_transcriptome_reads = matrix_to_use.sum()
pct_ebv_reads = (total_viral_reads / total_transcriptome_reads) * 100

print("\n" + "="*40)
print("   📊 EBV WT DATASET STATISTICS")
print("="*40)
print(f"Total EBV Genes found in Dataset: {len(ebv_genes)}")
print("-" * 40)
print(f"Total Cells:           {total_cells:,}")
print(f"Infected Cells (EBV+): {infected_cells:,}")
print(f"% EBV+ Cells:          {pct_ebv_cells:.4f}%")
print("-" * 40)
print(f"Total Viral Reads:     {int(total_viral_reads):,}")
print(f"Total Reads (All):     {int(total_transcriptome_reads):,}")
print(f"% of EBV Reads:        {pct_ebv_reads:.6f}%")
print("="*40 + "\n")
