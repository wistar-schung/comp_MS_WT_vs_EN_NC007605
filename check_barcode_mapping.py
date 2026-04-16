import scanpy as sc
import pandas as pd
import re
import os
import warnings

warnings.filterwarnings('ignore')

# ==========================================
# 1. Define File Paths
# ==========================================
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
OUTPUT_CSV = "Enriched_Barcode_WT_Mapping_Status.csv"

# ==========================================
# 2. Extract Metadata & Parse Barcodes
# ==========================================
print(f"Loading QC Data: {FILE_PATH} (in backed mode to save RAM)")
adata = sc.read_h5ad(FILE_PATH, backed='r') 

print("Extracting metadata and isolating core Parse Biosciences cell barcodes...")
df = pd.DataFrame(index=adata.obs_names)

# ---> SAFE PARSE BIOSCIENCES EXTRACTION <---
# Strips only the appended library identifiers from the end of the string, preserving the Parse combinatorial barcode

# ---> UPDATED SAFE EXTRACTION <---
# Now explicitly includes _ENR
df_obs_names = pd.Series(adata.obs_names)

# ---> PARSE BIOSCIENCES EXACT BARCODE MATCHING <---
# Ignore the integer row indices and pull the true combinatorial barcode directly
df['Core_Barcode'] = adata.obs['bc_wells'].astype(str).values

# Safely extract dataset and cell type
df['Dataset'] = adata.obs['dataset'].astype(str).values
df['Cell_Type'] = adata.obs['majority_voting'].astype(str).values

# Clean up conditions
df['Disease'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').astype(str).apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
)
df['Infection'] = adata.obs.get('Infection', 'Unknown').astype(str).apply(
    lambda x: 'Mock' if 'Mock' in str(x) else 'EBV'
)
df['Day_Sort'] = adata.obs.get('Day', 'Unknown').astype(str).apply(
    lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 999
)
df['Day'] = df['Day_Sort'].apply(lambda x: f"Day {x}" if x != 999 else "Unknown")

# ==========================================
# 3. Barcode-Level Intersection Logic
# ==========================================
print("Identifying WT vs Enriched mappings at the individual barcode level...")

# Create a unique key for the physical well + the barcode
df['Well_Barcode_Key'] = df['Day'].astype(str) + "_" + df['Disease'].astype(str) + "_" + df['Infection'].astype(str) + "_" + df['Core_Barcode'].astype(str)

# 1. Get the exact set of all valid WT Well+Barcode keys
wt_keys = set(df[df['Dataset'] == 'WT']['Well_Barcode_Key'])

# 2. Isolate just the Enriched cells
enr_df = df[df['Dataset'] == 'Enriched'].copy()

# 3. Safely map and build output
if enr_df.empty:
    print("\nWarning: No 'Enriched' cells found in the dataset.")
    output_df = pd.DataFrame()
else:
    enr_df['Mapped_to_WT'] = enr_df['Well_Barcode_Key'].isin(wt_keys)
    output_df = enr_df[['Day', 'Disease', 'Infection', 'Cell_Type', 'Core_Barcode', 'Mapped_to_WT']].copy()
    output_df['Original_adata_Index'] = enr_df.index

# ==========================================
# 4. Calculate Summary & Export
# ==========================================
total_enr = len(output_df) if not output_df.empty else 0
mapped_enr = output_df['Mapped_to_WT'].sum() if total_enr > 0 else 0
unmapped_enr = total_enr - mapped_enr

print("\n" + "="*80)
print(" 🔬 ENRICHED BARCODE MAPPING SUMMARY")
print("="*80)
if total_enr > 0:
    print(f"Total Enriched Cells: {total_enr:,}")
    print(f"Successfully Mapped to WT: {mapped_enr:,} ({(mapped_enr/total_enr*100):.2f}%)")
    print(f"Orphaned (Failed WT QC): {unmapped_enr:,} ({(unmapped_enr/total_enr*100):.2f}%)")
else:
    print("No Enriched cells to map.")
print("="*80 + "\n")

if total_enr > 0:
    output_df.to_csv(OUTPUT_CSV, index=False)
    print(f"✅ Detailed barcode-level mapping saved to: {os.path.abspath(OUTPUT_CSV)}")
