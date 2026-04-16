import scanpy as sc
import pandas as pd
import re
import os
import warnings

warnings.filterwarnings('ignore')

# ==========================================
# 1. Define File Paths
# ==========================================
FILE_PRE_QC = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_MULTI_compartment.h5ad"
FILE_POST_QC = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
OUTPUT_CSV = "QC_Dropout_Comparison_Report.csv"

# ==========================================
# 2. Helper Function: Clean Metadata
# ==========================================
def extract_clean_metadata(adata_obj, name):
    print(f"Extracting metadata from {name}...")
    df = pd.DataFrame(index=adata_obj.obs.index)
    
    # Safely extract dataset and force back to raw string to avoid Categorical errors
    df['Dataset'] = adata_obj.obs['dataset'].astype(str)
    
    # Clean Disease Group
    df['Disease'] = adata_obj.obs.get('Disease_Condition (Detail)', 'Unknown').astype(str).apply(
        lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
    )
    
    # Clean Infection
    df['Infection'] = adata_obj.obs.get('Infection', 'Unknown').astype(str).apply(
        lambda x: 'Mock' if 'Mock' in str(x) else 'EBV'
    )
    
    # Clean Day (extract integer and format nicely)
    def clean_day(day_str):
        match = re.search(r'\d+', str(day_str))
        return f"Day {match.group()}" if match else 'Unknown'
    
    df['Day'] = adata_obj.obs.get('Day', 'Unknown').astype(str).apply(clean_day)
    
    # Create sort key for chronological printing
    df['Day_Sort'] = adata_obj.obs.get('Day', 'Unknown').astype(str).apply(
        lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else 999
    )
    
    return df

# ==========================================
# 3. Load Data & Calculate Group Sizes
# ==========================================
print(f"Loading Pre-QC Data: {FILE_PRE_QC}")
adata_pre = sc.read_h5ad(FILE_PRE_QC, backed='r') 
df_pre = extract_clean_metadata(adata_pre, "Pre-QC")

print(f"Loading Post-QC Data: {FILE_POST_QC}")
adata_post = sc.read_h5ad(FILE_POST_QC, backed='r')
df_post = extract_clean_metadata(adata_post, "Post-QC")

# Group by the specific well conditions and count cells
group_cols = ['Day_Sort', 'Day', 'Disease', 'Infection', 'Dataset']

print("\nCalculating summary statistics...")
summary_pre = df_pre.groupby(group_cols).size().reset_index(name='Pre_QC_Cells')
summary_post = df_post.groupby(group_cols).size().reset_index(name='Post_QC_Cells')

# ==========================================
# 4. Merge & Calculate Dropout Metrics
# ==========================================
# Merge the two summaries together
comparison_df = pd.merge(summary_pre, summary_post, on=group_cols, how='outer')

# ---> BUG FIX: Fill NA ONLY on the specific integer columns <---
comparison_df['Pre_QC_Cells'] = comparison_df['Pre_QC_Cells'].fillna(0).astype(int)
comparison_df['Post_QC_Cells'] = comparison_df['Post_QC_Cells'].fillna(0).astype(int)

# Calculate Dropout Mathematics
comparison_df['Cells_Removed'] = comparison_df['Pre_QC_Cells'] - comparison_df['Post_QC_Cells']

# Safely calculate percentages (replacing 0 with 1 to avoid dividing by zero if a condition didn't exist)
comparison_df['%_Retained'] = (comparison_df['Post_QC_Cells'] / comparison_df['Pre_QC_Cells'].replace(0, 1) * 100).round(2)
comparison_df['%_Lost'] = (comparison_df['Cells_Removed'] / comparison_df['Pre_QC_Cells'].replace(0, 1) * 100).round(2)

# Sort logically for readability
comparison_df = comparison_df.sort_values(by=['Day_Sort', 'Disease', 'Infection', 'Dataset']).drop(columns=['Day_Sort'])

# ==========================================
# 5. Output Results
# ==========================================
print("\n" + "="*80)
print(" 🧬 QC DROPOUT COMPARISON REPORT: PRE-QC vs. POST-QC")
print("="*80)
print(comparison_df.to_string(index=False))
print("="*80 + "\n")

# Save to CSV
comparison_df.to_csv(OUTPUT_CSV, index=False)
print(f"✅ Detailed report saved to: {os.path.abspath(OUTPUT_CSV)}")
