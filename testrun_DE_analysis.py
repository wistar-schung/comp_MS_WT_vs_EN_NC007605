import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import os
import sys
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import scipy.sparse as sp


# ==========================================
# 0. CONFIGURATION & RUN MODES
# ==========================================
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"

# 🛑 STEP 1: Leave this as True first. Run the script. 
# It will print your metadata and immediately exit.
RUN_DISCOVERY_MODE = False 

# 🟢 STEP 2: Once you find your column name from the printout above, 
# set RUN_DISCOVERY_MODE to False, paste the name below, and run again.
BIOLOGICAL_REPLICATE_COL = "PatientID" # e.g., 'Donor', 'Patient_ID', 'bc_wells'

# Test Parameters for Execution Mode (Keep it small for a fast test)
TARGET_CELL_TYPE = "Memory B cells"
TEST_GROUP = "MS Active_EBV_Day 7"
REF_GROUP = "MS Active_EBV_Day 1"

# ==========================================
# 1. LOAD DATA 
# ==========================================
print(f"Loading data from {FILE_PATH}...")
# We load the full object, but we will subset it immediately depending on the mode
adata = sc.read_h5ad(FILE_PATH)

# ==========================================
# 2. DISCOVERY MODE (Find the Replicate Column)
# ==========================================
if RUN_DISCOVERY_MODE:
    print("\n" + "="*50)
    print("🕵️ METADATA DISCOVERY MODE ACTIVE")
    print("="*50)
    
    print("\n[1] ALL AVAILABLE METADATA COLUMNS:")
    all_cols = adata.obs.columns.tolist()
    print(all_cols)
    
    print("\n[2] HUNTING FOR LIKELY PATIENT/DONOR COLUMNS...")
    # Keywords that usually denote biological replicates
    keywords = ['donor', 'patient', 'sample', 'id', 'well', 'batch', 'rep']
    candidates = [col for col in all_cols if any(k in col.lower() for k in keywords)]
    
    if not candidates:
        print("Could not automatically guess the column. Please review the full list above.")
    else:
        for col in candidates:
            print(f"\n--- Inspecting Candidate Column: '{col}' ---")
            unique_count = adata.obs[col].nunique()
            print(f"Total Unique Values: {unique_count}")
            print(f"Top 5 most common values (and their cell counts):")
            print(adata.obs[col].value_counts().head(5).to_string())
            
            # If a column has thousands of unique values, it's a cell barcode, not a patient ID
            if unique_count > 500:
                print("⚠️ WARNING: High number of unique values. This is likely a cell barcode, NOT a patient identifier.")
    
    print("\n" + "="*50)
    print("ACTION REQUIRED:")
    print("1. Find the column name that uniquely identifies your original patients/donors.")
    print("2. Change 'RUN_DISCOVERY_MODE = False' at the top of this script.")
    print("3. Update 'BIOLOGICAL_REPLICATE_COL' with your chosen column name.")
    print("4. Re-run the script to execute the PyDESeq2 pipeline.")
    print("="*50 + "\n")
    sys.exit(0) # Stop the script here so we don't waste time running DE

# ==========================================
# 3. EXECUTION MODE (Pseudo-bulk & PyDESeq2)
# ==========================================
print("\n" + "="*50)
print(f"🚀 FAST DE EXECUTION MODE ACTIVE")
print(f"Using '{BIOLOGICAL_REPLICATE_COL}' as the biological replicate.")
print("="*50)

# Check if the column actually exists
if BIOLOGICAL_REPLICATE_COL not in adata.obs.columns:
    print(f"❌ ERROR: The column '{BIOLOGICAL_REPLICATE_COL}' does not exist in adata.obs.")
    sys.exit(1)

print(f"\nSubsetting to target cell type: {TARGET_CELL_TYPE}...")
adata_sub = adata[adata.obs['majority_voting'] == TARGET_CELL_TYPE].copy()


# --- NEW: SAFE RAW COUNT EXTRACTION ---
print("Securing raw integer counts for DESeq2...")

# 1. Recover from .raw if it exists
if adata_sub.raw is not None:
    print("  -> Recovering data from adata.raw")
    adata_sub = adata_sub.raw.to_adata()

# 2. If your pipeline saved raw counts in a layer, move it to .X
if 'counts' in adata_sub.layers:
    print("  -> Moving adata.layers['counts'] to .X")
    adata_sub.X = adata_sub.layers['counts'].copy()

# 3. Force integer casting. 
# DESeq2 strictly requires ints. We round first to fix precision errors (e.g., 1.0000001 -> 1)
print("  -> Casting matrix to strict integer dtype...")
if sp.issparse(adata_sub.X):
    # For sparse matrices
    adata_sub.X.data = np.round(adata_sub.X.data).astype(int)
else:
    # For dense matrices
    adata_sub.X = np.round(adata_sub.X).astype(int)
# --------------------------------------

# Create the specific DE Interaction Column
print("Generating DE_Combo (Disease_Infection_Day) column...")
adata_sub.obs['Disease_Group'] = adata_sub.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else 'HC')
)
adata_sub.obs['Infection_Clean'] = adata_sub.obs.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV')
adata_sub.obs['Day_Clean'] = adata_sub.obs['Day'].astype(str)

adata_sub.obs['DE_Combo'] = (
    adata_sub.obs['Disease_Group'].astype(str) + "_" + 
    adata_sub.obs['Infection_Clean'].astype(str) + "_" + 
    adata_sub.obs['Day_Clean'].astype(str)
)

# Crucial Speed Step: Subset down to ONLY the two groups we are comparing right now
print(f"Filtering to only keep: '{TEST_GROUP}' and '{REF_GROUP}'...")
adata_test = adata_sub[adata_sub.obs['DE_Combo'].isin([TEST_GROUP, REF_GROUP])].copy()

if adata_test.n_obs < 10:
    print("❌ ERROR: Not enough cells found for these test groups. Check your naming conventions.")
    sys.exit(1)



# Force PatientID to be a string so decoupler can concatenate it
adata_test.obs[BIOLOGICAL_REPLICATE_COL] = adata_test.obs[BIOLOGICAL_REPLICATE_COL].astype(str)

# Generate Pseudo-bulk using the new decoupler 2.0+ API
print(f"Pseudo-bulking by {BIOLOGICAL_REPLICATE_COL}...")
pbulk = dc.pp.pseudobulk(
    adata_test,
    sample_col=BIOLOGICAL_REPLICATE_COL,
    groups_col='DE_Combo',
    mode='sum'
)

# Optional quality control: Filter out patient/condition combos that didn't have enough cells
# dc.pp.filter_samples(pbulk, min_cells=10, min_counts=1000)

# Extract counts and ensure it is a dense Pandas DataFrame for PyDESeq2
count_data = pbulk.X.toarray() if sp.issparse(pbulk.X) else pbulk.X
counts = pd.DataFrame(count_data, index=pbulk.obs_names, columns=pbulk.var_names)

# Extract metadata
metadata = pbulk.obs[['DE_Combo', BIOLOGICAL_REPLICATE_COL]].copy()

# Verify Replicate Math
test_n = sum(metadata['DE_Combo'] == TEST_GROUP)
ref_n = sum(metadata['DE_Combo'] == REF_GROUP)
print(f"Replicates found -> {TEST_GROUP}: n={test_n} | {REF_GROUP}: n={ref_n}")

if test_n < 2 or ref_n < 2:
    print("\n❌ CRITICAL FAILURE: DESeq2 requires at least 2 independent donors per group.")
    print("Pseudo-bulking cannot be performed on this specific contrast because there are not enough biological replicates.")
    sys.exit(1)

# Initialize and Run PyDESeq2
print("\nInitializing PyDESeq2...")
dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design_factors='DE_Combo'
)

print("Fitting DESeq2 model (Dispersion estimation)...")
dds.deseq2()

print("Extracting Stats...")
stat_res = DeseqStats(
    dds,
    contrast=['DE_Combo', TEST_GROUP, REF_GROUP]
)
stat_res.summary()

# Format and Save Results
res_df = stat_res.results_df.dropna(subset=['padj']).reset_index()
res_df.rename(columns={'index': 'Gene'}, inplace=True)
res_df = res_df.sort_values('padj')

output_file = f"Test_DE_{TARGET_CELL_TYPE.replace(' ', '_')}_{TEST_GROUP}_vs_{REF_GROUP}.csv".replace(" ", "")
res_df.to_csv(output_file, index=False)

print(f"\n✅ SUCCESS! Fast test complete.")
print(f"Top 5 DE Genes:")
print(res_df[['Gene', 'log2FoldChange', 'padj']].head(5).to_string(index=False))
print(f"\nFull results saved to: {output_file}")
