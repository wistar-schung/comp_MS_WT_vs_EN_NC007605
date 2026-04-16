import scanpy as sc
import pandas as pd
import os
import re

FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"

print("Loading QC h5ad dataset...")
adata = sc.read_h5ad(FILE_PATH)

print("Loading Metadata.csv...")
meta = pd.read_csv("Metadata.csv")
meta.columns = meta.columns.str.strip()

# ==========================================
# 1. Merge the Disease Columns
# ==========================================
# Combines "MS" and "Active" into "MS Active", handles HC's "NA" gracefully
def get_full_disease(row):
    disease = str(row['Disease']).strip()
    cond = str(row['Disease_condition']).strip()
    if disease == 'HC': 
        return 'HC'
    if disease == 'MS':
        if cond == 'Active': return 'MS Active'
        if cond == 'Stable': return 'MS Stable'
    return 'Unknown'
    
meta['Full_Disease'] = meta.apply(get_full_disease, axis=1)

# ==========================================
# 2. Format the Day Strings
# ==========================================
# Converts "D1" to "Day 1" so the DE script can find the cells
def clean_day(d):
    m = re.search(r'\d+', str(d))
    return f"Day {m.group()}" if m else "Unknown"
    
meta['Clean_Day'] = meta['Day'].apply(clean_day)

# ==========================================
# 3. Apply the Mappings
# ==========================================
well_to_infection = dict(zip(meta['Well'], meta['Infection']))
well_to_day = dict(zip(meta['Well'], meta['Clean_Day']))
well_to_disease = dict(zip(meta['Well'], meta['Full_Disease']))

print("Applying strict physical well mappings to cell barcodes...")
adata.obs['Infection'] = adata.obs['bc1_well'].map(well_to_infection).fillna('Unknown')
adata.obs['Day'] = adata.obs['bc1_well'].map(well_to_day).fillna('Unknown')
adata.obs['Disease_Condition (Detail)'] = adata.obs['bc1_well'].map(well_to_disease).fillna('Unknown')

print("Saving corrected dataset...")
adata.write_h5ad(FILE_PATH)
print("✅ Metadata perfectly aligned! You can now re-run generate_master_dashboard.py")