import scanpy as sc
import pandas as pd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest

# 1. Setup Data
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
print(f"[{time.strftime('%H:%M:%S')}] Loading {FILE_PATH}...")
start_time = time.time()
adata = sc.read_h5ad(FILE_PATH)
print(f"[{time.strftime('%H:%M:%S')}] Data loaded in {time.time() - start_time:.2f}s. Shape: {adata.shape}")

# Clean up metadata
print(f"[{time.strftime('%H:%M:%S')}] Cleaning metadata...")
adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
adata.var_names_make_unique()
adata.obs['Disease_Group'] = adata.obs['Disease_Condition (Detail)'].astype(str)
adata.obs['Infection_Clean'] = adata.obs['Infection'].astype(str)
adata.obs['Day_Clean'] = adata.obs['Day'].astype(str)
adata.obs['Analysis_Group'] = adata.obs['Disease_Group'] + "_" + adata.obs['Day_Clean'] + "_" + adata.obs['Infection_Clean']

# ==============================================================================
# PART 1: Marker Gene Overlap Analysis
# ==============================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] --- Identifying Marker Genes per Cell Type ---")
marker_start = time.time()
sc.tl.rank_genes_groups(adata, groupby='majority_voting', method='t-test', n_genes=100)
print(f"[{time.strftime('%H:%M:%S')}] Marker identification complete in {time.time() - marker_start:.2f}s.")

cell_types = sorted(adata.obs['majority_voting'].unique().tolist())
marker_dict_top = {ct: set(adata.uns['rank_genes_groups']['names'][ct][:10]) for ct in cell_types}

print(f"[{time.strftime('%H:%M:%S')}] --- Generating Marker Overlap Heatmap (Threshold > 0.5) ---")
overlap_matrix = np.zeros((len(cell_types), len(cell_types)))
overlaps = []

for i, ct1 in enumerate(cell_types):
    for j, ct2 in enumerate(cell_types):
        shared = marker_dict_top[ct1].intersection(marker_dict_top[ct2])
        ratio = len(shared) / 10.0
        overlap_matrix[i, j] = ratio
        if i < j and ratio >= 0.5:
            overlaps.append({'CT1': ct1, 'CT2': ct2, 'Overlap_Ratio': ratio, 'Genes': list(shared)})

# Visualizing Overlap
plt.figure(figsize=(12, 10))
sns.heatmap(overlap_matrix, xticklabels=cell_types, yticklabels=cell_types, annot=True, fmt=".1f", cmap="YlOrRd")
plt.title("Marker Overlap Ratio (Top 10 Genes)")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig("marker_overlap_heatmap.png")
plt.close()

overlap_df = pd.DataFrame(overlaps)
overlap_df.to_csv("cell_type_marker_overlaps_v2.csv", index=False)

# ==============================================================================
# PART 2: Detailed Population Composition Analysis
# ==============================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] --- Analyzing Population Changes across Groups ---")

# Calculate proportions (using WT only for cleaner stats)
wt_mask = adata.obs['dataset'] == 'WT'
pop_counts = adata.obs[wt_mask].groupby(['Analysis_Group', 'majority_voting']).size().unstack(fill_value=0)
group_totals = pop_counts.sum(axis=1)
pop_pct = pop_counts.div(group_totals, axis=0) * 100

# 2.1 Disease Group Variation (HC vs MS Active vs MS Stable)
baseline_groups = [g for g in pop_pct.index if "Day 1" in g and "Mock" in g]
disease_baseline = pop_pct.loc[baseline_groups].copy()
disease_baseline.index = [i.split("_")[0] for i in disease_baseline.index]

print(f"[{time.strftime('%H:%M:%S')}] Identifying cell types most different among disease groups...")
disease_diff = disease_baseline.std().sort_values(ascending=False).to_frame(name='Std_Dev_Across_Diseases')
disease_diff['Max_Diff'] = disease_baseline.max() - disease_baseline.min()
disease_diff.to_csv("disease_group_composition_variation.csv")

# 2.2 Temporal Changes & Range Filtering
print(f"[{time.strftime('%H:%M:%S')}] Analyzing temporal changes (Filtering Range > 10)...")
temporal_res = []
for disease in ['HC', 'MS Active', 'MS Stable']:
    for infection in ['Mock', 'EBV']:
        sub_groups = [g for g in pop_pct.index if g.startswith(disease) and g.endswith(infection)]
        if len(sub_groups) < 2: continue
        
        sub_pct = pop_pct.loc[sub_groups]
        for ct in cell_types:
            r = sub_pct[ct].max() - sub_pct[ct].min()
            temporal_res.append({
                'Disease': disease, 'Infection': infection, 'Cell_Type': ct, 'Range': r
            })

temporal_df = pd.DataFrame(temporal_res)
high_range_cts = temporal_df[temporal_df['Range'] > 10]['Cell_Type'].unique()

# Identify "Double Dynamic" Cell Types (Dynamic in BOTH Mock and EBV for same disease)
double_dynamic = []
for disease in ['HC', 'MS Active', 'MS Stable']:
    mock_high = set(temporal_df[(temporal_df['Disease'] == disease) & (temporal_df['Infection'] == 'Mock') & (temporal_df['Range'] > 10)]['Cell_Type'])
    ebv_high = set(temporal_df[(temporal_df['Disease'] == disease) & (temporal_df['Infection'] == 'EBV') & (temporal_df['Range'] > 10)]['Cell_Type'])
    shared = mock_high.intersection(ebv_high)
    for ct in shared:
        double_dynamic.append({'Disease': disease, 'Cell_Type': ct, 'Status': 'Highly Dynamic in BOTH Mock & EBV'})

pd.DataFrame(double_dynamic).to_csv("double_dynamic_cell_types.csv", index=False)

# 2.3 Relative Mock-Normalized EBV Response: (EBV - Mock) / Mock
print(f"[{time.strftime('%H:%M:%S')}] Calculating Relative (Fold Change) EBV response...")
ebv_shifts = []
epsilon = 1e-6 # Prevent division by zero
for disease in ['HC', 'MS Active', 'MS Stable']:
    for day in ['Day 1', 'Day 7', 'Day 15']:
        mock_g = f"{disease}_{day}_Mock"; ebv_g = f"{disease}_{day}_EBV"
        if mock_g in pop_pct.index and ebv_g in pop_pct.index:
            # (EBV - Mock) / Mock
            rel_shift = (pop_pct.loc[ebv_g] - pop_pct.loc[mock_g]) / (pop_pct.loc[mock_g] + epsilon)
            for ct, val in rel_shift.items():
                ebv_shifts.append({'Disease': disease, 'Day': day, 'Cell_Type': ct, 'Relative_Shift': val})

shift_df = pd.DataFrame(ebv_shifts)
shift_df.to_csv("ebv_relative_composition_shifts.csv", index=False)

# ==============================================================================
# PART 3: Visualizations
# ==============================================================================
print(f"\n[{time.strftime('%H:%M:%S')}] --- Generating Composition Visualizations ---")

# 3.1 MASTER HEATMAP: Z-Score per Cell Type
print(f"[{time.strftime('%H:%M:%S')}] Generating Z-score Heatmap...")
pop_pct_z = pop_pct.apply(lambda x: (x - x.mean()) / (x.std() + epsilon), axis=0)
plt.figure(figsize=(18, 12))
sns.heatmap(pop_pct_z.T, cmap="RdBu_r", center=0, annot=False, cbar_kws={'label': 'Z-Score (Relative Abundance)'})
plt.title("Z-Score Normalized Cell Composition (Relative Changes within each Cell Type)")
plt.tight_layout()
plt.savefig("global_composition_zscore_heatmap.png")
plt.close()

# 3.2 High-Dynamics Visualization (Range > 10)
if len(high_range_cts) > 0:
    print(f"[{time.strftime('%H:%M:%S')}] Plotting High-Dynamics Cell Types...")
    dyn_plot_df = temporal_df[temporal_df['Cell_Type'].isin(high_range_cts)]
    plt.figure(figsize=(14, 8))
    sns.barplot(data=dyn_plot_df, x='Cell_Type', y='Range', hue='Disease')
    plt.axhline(10, color='red', linestyle='--', label='Threshold (Range=10)')
    plt.title("Cell Types with High Temporal Dynamics (Range > 10%)")
    plt.xticks(rotation=45, ha='right')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("high_range_temporal_dynamics.png")
    plt.close()

# 3.3 Relative EBV Shift Heatmap
pivot_shift_rel = shift_df.pivot_table(index='Cell_Type', columns=['Disease', 'Day'], values='Relative_Shift')
plt.figure(figsize=(16, 12))
sns.heatmap(pivot_shift_rel, cmap="PRGn", center=0, annot=False)
plt.title("Relative EBV Response: (EBV - Mock) / Mock")
plt.tight_layout()
plt.savefig("ebv_relative_shift_heatmap.png")
plt.close()

print(f"\n[{time.strftime('%H:%M:%S')}] Analysis complete. Results saved to:")
print("- global_composition_zscore_heatmap.png (Small populations visible)")
print("- high_range_temporal_dynamics.png (Range > 10 filtered)")
print("- ebv_relative_shift_heatmap.png (Fold-change normalization)")
print("- double_dynamic_cell_types.csv (Highlighting Mock & EBV shifts)")

print(f"\n[{time.strftime('%H:%M:%S')}] Analysis complete. Results saved to:")
print("- marker_overlap_heatmap.png")
print("- global_composition_heatmap.png")
print("- top_disease_differences.png")
print("- ebv_shift_heatmap.png")
print("- disease_group_composition_variation.csv")
print("- temporal_composition_changes.csv")
print("- ebv_induced_composition_shifts.csv")
print(f"Total time: {time.time() - start_time:.2f}s.")
