import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.sparse as sp
import os
import re
import warnings

warnings.filterwarnings('ignore')

# ==========================================
# 1. Setup & Configuration
# ==========================================
print("Loading combined QC dataset...")
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
adata = sc.read_h5ad(FILE_PATH)

OUTPUT_DIR = "Publication_UMAPs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Clean var names 
adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
adata.var_names_make_unique()

# ==========================================
# 2. Strict Viral Load Calculation (GTF)
# ==========================================
print("Calculating exact viral loads from raw counts...")
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

ebv_genes = [g for g in adata.var_names if g in EXACT_EBV_GENES or g.startswith('EBV-')]

if adata.raw is not None:
    raw_var_names = adata.raw.var_names.str.replace(r'(?i)_type1$', '', regex=True)
    ebv_raw_indices = [i for i, g in enumerate(raw_var_names) if g in ebv_genes]
    ebv_matrix = adata.raw.X[:, ebv_raw_indices]
else:
    ebv_matrix = adata[:, ebv_genes].X

if sp.issparse(ebv_matrix):
    adata.obs['Viral_Counts'] = np.ravel(ebv_matrix.sum(axis=1))
else:
    adata.obs['Viral_Counts'] = np.ravel(np.sum(ebv_matrix, axis=1))

# ==========================================
# 3. Metadata & Barcode Standardization
# ==========================================
CELLTYPE_COLOR_MAP = {
    'Tcm/Naive helper T cells': '#1f77b4', 'Tem cytotoxic T cells': '#2ca02c',
    'Plasmablasts': '#d62728', 'NK cells': '#9467bd', 
    'Regulatory T cells': '#98df8a', 'Classical monocytes': '#ff9896',
    'Memory B cells': '#fdbf6f', 'Macrophages': '#c5b0d5',
    'Tem/Effector helper T cells PD1+': '#f47d4d', 'Cycling T cells': '#75a3e1',
    'Age-associated B cells': '#299d5c', 'Naive B cells': '#9e9bc1',
    'CD16- NK cells': '#f1a243', 'Migratory DCs': '#7196aa',
    'Plasma cells': '#8bcebb', 'Innate Lymphoid Cells': '#1e9a79', 
    'MAIT cells': '#afd85c', 'CD16+ NK cells': '#fadd4d',
    'Non-classical monocytes': '#c95b1a', 'pDC': '#786db1', 'Unknown': '#e3e3e3'
}

def get_day_sort_key(day_str):
    match = re.search(r'\d+', str(day_str))
    return int(match.group()) if match else 9999

adata.obs['Day_Sort'] = adata.obs.get('Day', 'Unknown').astype(str).apply(get_day_sort_key)
adata.obs['Disease_Group'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
)
adata.obs['Infection_Clean'] = adata.obs.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV')

# ---> PARSE BIOSCIENCES EXACT BARCODE MATCHING <---
# Use the true split-pool well combinations to link WT to Enriched cells
adata.obs['Core_Barcode'] = adata.obs['bc_wells'].astype(str).values

TARGET_DAYS = {1: 'Day 1', 7: 'Day 7', 15: 'Day 15'}
DISEASES = ['HC', 'MS Active', 'MS Stable']
INFECTIONS = ['Mock', 'EBV']

umap_coords = adata.obsm['X_umap']
global_wt_mask = adata.obs['dataset'] == 'WT'

# ==========================================
# 4. Generate Intersected Publication Plots
# ==========================================
print("\nGenerating High-Resolution Condition Overlays with Strict Barcode Intersection...")

summary_stats = []

for day_int, day_str in TARGET_DAYS.items():
    for disease in DISEASES:
        for inf in INFECTIONS:
            
            cond_name = f"{day_str} {disease} {inf}"
            safe_name = cond_name.replace(" ", "_")
            
            # 1. Isolate the base WT cells for this exact condition
            mask_cond_wt = (adata.obs['dataset'] == 'WT') & (adata.obs['Day_Sort'] == day_int) & (adata.obs['Disease_Group'] == disease) & (adata.obs['Infection_Clean'] == inf)
            
            # 2. Extract the core barcodes of the WT cells that survived QC
            wt_core_barcodes = set(adata.obs.loc[mask_cond_wt, 'Core_Barcode'])
            
            # 3. Find the Enriched cells for this condition with UMI >= 2...
            mask_cond_enr_ebv_raw = (adata.obs['dataset'] == 'Enriched') & (adata.obs['Day_Sort'] == day_int) & (adata.obs['Disease_Group'] == disease) & (adata.obs['Infection_Clean'] == inf) & (adata.obs['Viral_Counts'] >= 2)
            
            # 4. STRICT INTERSECTION: Only keep the Enriched cell if its twin WT barcode exists in the base layer!
            mask_cond_enr_ebv_intersected = mask_cond_enr_ebv_raw & adata.obs['Core_Barcode'].isin(wt_core_barcodes)
            
            num_wt = mask_cond_wt.sum()
            num_ebv_raw = mask_cond_enr_ebv_raw.sum()
            num_ebv_final = mask_cond_enr_ebv_intersected.sum()
            
            # Save to summary table
            summary_stats.append({
                'Day': day_str, 'Disease': disease, 'Infection': inf,
                'WT_Base_Cells_Plotted': num_wt,
                'Enriched_EBV_Cells_Raw': num_ebv_raw,
                'Enriched_EBV_Cells_Plotted_Intersected': num_ebv_final
            })
            
            # Skip if there is virtually no WT base data for this well
            if num_wt < 20:
                continue
                
            print(f" -> Plotting {cond_name} (Base: {num_wt} WT cells | Intersected EBV+ Dots: {num_ebv_final})")
            
            fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
            
            # Layer A: Universal Background (All WT cells globally, light grey, to orient the global shape)
            ax.scatter(umap_coords[global_wt_mask, 0], umap_coords[global_wt_mask, 1], c='#E8E8E8', s=1, alpha=0.3, edgecolors='none', zorder=1)
            
            # Layer B: Condition-Specific Base (WT cells colored by exact Cell Type)
            cond_wt_indices = np.where(mask_cond_wt)[0]
            cell_types = adata.obs['majority_voting'].iloc[cond_wt_indices].astype(str)
            colors = [CELLTYPE_COLOR_MAP.get(ct, "#d3d3d3") for ct in cell_types]
            
            ax.scatter(umap_coords[cond_wt_indices, 0], umap_coords[cond_wt_indices, 1], 
                       c=colors, s=5, alpha=0.9, edgecolors='none', zorder=2)
                       
            # Layer C: Intersected Enriched EBV+ Cells (Bright Ruby Red, with dark borders)
            if num_ebv_final > 0:
                cond_ebv_indices = np.where(mask_cond_enr_ebv_intersected)[0]
                ax.scatter(umap_coords[cond_ebv_indices, 0], umap_coords[cond_ebv_indices, 1], 
                           c='#e6194b', s=25, alpha=1.0, edgecolors='#4a0000', linewidths=0.5, zorder=3)

            # --- Dynamic Legend Generation ---
            legend_handles = []
            unique_cts = adata.obs['majority_voting'].iloc[cond_wt_indices].unique()
            for ct in sorted(unique_cts):
                if ct in CELLTYPE_COLOR_MAP and str(ct) != 'Unknown':
                    legend_handles.append(mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor=CELLTYPE_COLOR_MAP[ct], markersize=8, label=ct))
            
            if num_ebv_final > 0:
                legend_handles.append(mlines.Line2D([0], [0], marker='o', color='w', markerfacecolor='#e6194b', markeredgecolor='#4a0000', markersize=10, label='EBV+ (UMI ≥ 2)'))

            if legend_handles:
                ax.legend(handles=legend_handles, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=11, title="Populations", title_fontsize=13)

            # --- Add Exact Counts to the UMAP Plot ---
            count_text = f"Total WT Cells: {num_wt:,}\nOverlapping EBV+ (UMI ≥ 2): {num_ebv_final:,}"
            ax.text(0.02, 0.98, count_text, transform=ax.transAxes, fontsize=12, fontweight='bold',
                    verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='#cccccc'))

            ax.set_title(f"{cond_name}", fontsize=18, fontweight='bold', pad=15)
            ax.axis('off') 
            
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, f"UMAP_{safe_name}.png"), dpi=300, bbox_inches='tight', transparent=False)
            plt.savefig(os.path.join(OUTPUT_DIR, f"UMAP_{safe_name}.pdf"), dpi=300, bbox_inches='tight', transparent=True)
            plt.close(fig)

# Save the final exact counts to a CSV table
summary_df = pd.DataFrame(summary_stats)
summary_path = os.path.join(OUTPUT_DIR, "UMAP_Plotted_Cell_Counts_Summary.csv")
summary_df.to_csv(summary_path, index=False)

print(f"\n✅ All plots saved to: {os.path.abspath(OUTPUT_DIR)}")
print(f"📊 Summary Table exported to: {summary_path}")
