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
OUTPUT_CSV = "WT_vs_Enriched_Barcode_Overlap_by_CellType.csv"

# ==========================================
# 2. Extract Metadata & Core Barcodes
# ==========================================
print(f"Loading QC Data: {FILE_PATH} (in backed mode to save RAM)")
adata = sc.read_h5ad(FILE_PATH, backed='r') 

print("Extracting metadata and isolating core 10X cell barcodes...")
df = pd.DataFrame(index=adata.obs_names)

# Extract core barcode to match physical droplets
df['Core_Barcode'] = pd.Series(adata.obs_names).str.split('-').str[0].str.split('_').str[0].values
df['Dataset'] = adata.obs['dataset'].values
df['Cell_Type'] = adata.obs['majority_voting'].values

# Clean up conditions
df['Disease'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
)
df['Infection'] = adata.obs.get('Infection', 'Unknown').apply(
    lambda x: 'Mock' if 'Mock' in str(x) else 'EBV'
)
df['Day_Sort'] = adata.obs.get('Day', 'Unknown').astype(str).apply(
    lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 999
)
df['Day'] = df['Day_Sort'].apply(lambda x: f"Day {x}" if x != 999 else "Unknown")

# ==========================================
# 3. Fast Grouping & Intersection Logic
# ==========================================
print("Calculating exact overlaps between WT and Enriched libraries by Cell Type...")

# Group by all experimental variables PLUS the exact cell type
grouping_cols = ['Day_Sort', 'Day', 'Disease', 'Infection', 'Cell_Type']

# Fast extraction of sets per condition
grouped = df.groupby(grouping_cols + ['Dataset'])['Core_Barcode'].apply(set).reset_index()

# Pivot the data so WT and Enriched sets are on the same row for easy intersection
pivot_df = grouped.pivot(index=grouping_cols, columns='Dataset', values='Core_Barcode').reset_index()

# Handle missing datasets gracefully (e.g., if a cell type has 0 WT cells)
if 'WT' not in pivot_df.columns: pivot_df['WT'] = [set() for _ in range(len(pivot_df))]
if 'Enriched' not in pivot_df.columns: pivot_df['Enriched'] = [set() for _ in range(len(pivot_df))]

# Fill NaNs with empty sets
pivot_df['WT'] = pivot_df['WT'].apply(lambda x: x if isinstance(x, set) else set())
pivot_df['Enriched'] = pivot_df['Enriched'].apply(lambda x: x if isinstance(x, set) else set())

# ==========================================
# 4. Calculate Mathematical Overlaps
# ==========================================
results = []

for _, row in pivot_df.iterrows():
    wt_barcodes = row['WT']
    enr_barcodes = row['Enriched']
    
    overlap = wt_barcodes.intersection(enr_barcodes)
    
    num_wt = len(wt_barcodes)
    num_enr = len(enr_barcodes)
    num_overlap = len(overlap)
    
    # Skip entirely empty groups
    if num_wt == 0 and num_enr == 0:
        continue
        
    pct_enr_in_wt = (num_overlap / num_enr * 100) if num_enr > 0 else 0
    pct_wt_in_enr = (num_overlap / num_wt * 100) if num_wt > 0 else 0
    
    results.append({
        'Day_Sort': row['Day_Sort'],
        'Day': row['Day'],
        'Disease': row['Disease'],
        'Infection': row['Infection'],
        'Cell_Type': row['Cell_Type'],
        'Total_WT_Cells': num_wt,
        'Total_Enriched_Cells': num_enr,
        'Overlapping_Cells': num_overlap,
        '%_of_Enriched_Found_in_WT': round(pct_enr_in_wt, 2),
        '%_of_WT_Found_in_Enriched': round(pct_wt_in_enr, 2)
    })

results_df = pd.DataFrame(results).sort_values(['Day_Sort', 'Disease', 'Infection', 'Cell_Type']).drop(columns=['Day_Sort'])

# ==========================================
# 5. Print & Export Results
# ==========================================
print("\n" + "="*120)
print(" 🔗 WT vs ENRICHED OVERLAPS BY CELL TYPE (POST-QC)")
print("="*120)
# Print the first 50 rows so the console doesn't get utterly spammed, but export everything
print(results_df.head(50).to_string(index=False))
print("... (Showing first 50 rows. Full data exported to CSV)")
print("="*120 + "\n")

results_df.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Detailed granular report saved to: {os.path.abspath(OUTPUT_CSV)}")
