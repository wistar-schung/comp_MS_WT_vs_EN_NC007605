import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import os
import mygene
import warnings

warnings.filterwarnings('ignore')

# ==========================================
# 0. Setup & Parameters
# ==========================================
sc.settings.verbosity = 3  
sc.settings.set_figure_params(dpi=150, facecolor='white')

FILE_WT = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/optimized_data/data_200_annotated_post_qc.h5ad"
FILE_ENR = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/optimized_data/data_EBV_enriched_100_post_qc.h5ad"

CELL_TYPE_COL = 'majority_voting' 

# Broad CellTypist Classes (Updated for PBMC biology)
B_CELL_CLASSES = ['Memory B cells', 'Plasma cells', 'Naive B cells', 'Age-associated B cells', 'Plasmablasts']
T_CELL_CLASSES = ['Tem cytotoxic T cells', 'Tem/Effector helper T cells PD1+', 'Tcm/Naive helper T cells', 'Regulatory T cells', 'Cycling T cells']
MONO_MAC_CLASSES = ['Classical monocytes', 'Non-classical monocytes', 'Macrophages']

TARGET_CELL_TYPES = B_CELL_CLASSES + T_CELL_CLASSES + MONO_MAC_CLASSES

# Targeted Gene Modules
MODULE_GENES = {
    'ABC_Score': ["ITGAX", "TBX21", "FCRL5", "ZEB2", "CXCR3", "LILRB1", "LILRB2"],
    'T_Exhaustion_Score': ["PDCD1", "LAG3", "HAVCR2", "CTLA4", "TIGIT", "TOX", "GZMB", "PRF1", "IFNG"],
    'Mono_Inflammatory_Score': ["IL1B", "TNF", "CXCL8", "CCL2", "NLRP3", "TLR4", "S100A8", "S100A9"],
    'IFN_Score': ["ISG15", "IFI44L", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "IRF7", "STAT1", "BST2", "XAF1"],
    'ROS_Curated_Score': ["SOD1", "SOD2", "SOD3", "CAT", "PRDX1", "PRDX2", "PRDX6", "GPX1", "GPX4", "GSR", "TXN", "TXNRD1", "NFE2L2", "HMOX1", "NQO1"]
}

KNOWN_EBV_GENES = [
    "A73", "BALF1", "BALF2", "BALF3", "BALF4", "BALF5", "BARF1", "BBLF1", "BBLF2/BBLF3", "BBLF4",
    "BBRF1", "BBRF2", "BBRF3", "BCRF1", "BDLF1", "BDLF2", "BDLF3", "BDLF3.5", "BDLF4", "BFLF1",
    "BFLF2", "BFRF1", "BFRF1A", "BFRF2", "BFRF3", "BGLF1", "BGLF2", "BGLF3", "BGLF3.5", "BGLF4",
    "BGLF5", "BGRF1/BDRF1", "BHRF1", "BILF1", "BILF2", "BKRF2", "BKRF3", "BKRF4", "BLLF1", "BLLF2",
    "BLLF3", "BLRF1", "BLRF2", "BMRF1", "BMRF2", "BNLF2a", "BNLF2b", "BNRF1", "BORF1", "BORF2",
    "BPLF1", "BRLF1", "BRRF1", "BRRF2", "BSLF1", "BSLF2/BMLF1", "BSRF1", "BTRF1", "BVLF1", "BVRF1",
    "BVRF2", "BXLF1", "BXLF2", "BXRF1", "BZLF1", "BZLF2", "BaRF1", "BcLF1", "BcRF1", "BdRF1",
    "EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B/EBNA-3C", "EBNA-LP", "LF1", "LF2", "LMP-1", "LMP-2A",
    "LMP-2B", "RPMS1"
]

DOWNSAMPLE = False         
MAX_CELLS = 1000000         

# ==========================================
# 1. Loading & Fixing Gene IDs
# ==========================================
print("\n--- 1. Loading Datasets ---")
adata_wt_full = sc.read_h5ad(FILE_WT, backed='r')
adata_enr = sc.read_h5ad(FILE_ENR)

# Strip hidden barcode suffixes to ensure perfect mapping
adata_wt_full.obs_names = adata_wt_full.obs_names.str.replace(r'-1$', '', regex=True)
adata_enr.obs_names = adata_enr.obs_names.str.replace(r'-1$', '', regex=True)

# --- Translate Enriched ENSG to Symbols ---
ensg_genes = [g for g in adata_enr.var_names if g.startswith('ENSG')]
if len(ensg_genes) > 0:
    print(f"Translating {len(ensg_genes)} ENSG IDs in Enriched data to match WT Gene Symbols...")
    mg = mygene.MyGeneInfo()
    clean_ids = [g.split('.')[0] for g in ensg_genes]
    
    query_results = mg.querymany(clean_ids, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True, verbose=False)
    
    ensembl_to_symbol = {}
    for orig, clean in zip(ensg_genes, clean_ids):
        if clean in query_results.index and 'symbol' in query_results.columns:
            sym = query_results.loc[clean, 'symbol']
            if isinstance(sym, pd.Series): sym = sym.iloc[0]
            if pd.notna(sym): ensembl_to_symbol[orig] = str(sym).upper()
            
    # Apply mapping, keeping viral genes untouched
    adata_enr.var_names = [ensembl_to_symbol.get(g, g) for g in adata_enr.var_names]
    adata_enr.var_names_make_unique()
    print("Translation complete.")

# --- Relabel Tissue-Resident/Non-PBMC Annotations ---
print("\nRe-labeling non-circulating immune populations to PBMC-relevant categories...")
relabel_map = {
    'Alveolar macrophages': 'Alveolar macrophages',
    'Intestinal macrophages': 'Intestinal macrophages',
    'ILC3': 'Innate Lymphoid Cells',
    'Tem/Trm cytotoxic T cells': 'Tem cytotoxic T cells'
}

if CELL_TYPE_COL in adata_wt_full.obs.columns:
    adata_wt_full.obs[CELL_TYPE_COL] = adata_wt_full.obs[CELL_TYPE_COL].replace(relabel_map)
if CELL_TYPE_COL in adata_enr.obs.columns:
    adata_enr.obs[CELL_TYPE_COL] = adata_enr.obs[CELL_TYPE_COL].replace(relabel_map)

# ==========================================
# 2. Label Transfer & Viral Bridging
# ==========================================
print("\n--- 2. Transferring WT Cell Annotations & Bridging Viral Data ---")
wt_label_dict = adata_wt_full.obs[CELL_TYPE_COL].to_dict()
adata_enr.obs['Original_Enriched_Prediction'] = adata_enr.obs.get(CELL_TYPE_COL, 'Unknown')

adata_enr.obs[CELL_TYPE_COL] = adata_enr.obs_names.map(
    lambda x: wt_label_dict.get(x, adata_enr.obs.loc[x, 'Original_Enriched_Prediction'])
)
print("Successfully anchored Enriched cells to their WT identities.")

# Robust Viral Gene Detection
ebv_genes = [g for g in adata_enr.var_names if g.split('_')[0] in KNOWN_EBV_GENES or 'TYPE1' in g.upper() or 'EBV' in g.upper()]

if ebv_genes:
    adata_enr.obs['Enriched_EBV_Counts'] = adata_enr[:, ebv_genes].X.sum(axis=1)
else:
    print("WARNING: No EBV genes found in Enriched dataset.")
    adata_enr.obs['Enriched_EBV_Counts'] = 0

enr_ebv_dict = adata_enr.obs['Enriched_EBV_Counts'].to_dict()
adata_wt_full.obs['Mapped_Enriched_EBV_Counts'] = adata_wt_full.obs_names.map(enr_ebv_dict).fillna(0)
print(f"Bridged viral data. Found {sum(adata_wt_full.obs['Mapped_Enriched_EBV_Counts'] > 0)} WT cells with confirmed Enriched viral loads.")

if DOWNSAMPLE and adata_wt_full.n_obs > MAX_CELLS:
    print(f"Extracting {MAX_CELLS:,} cells straight from disk into RAM...")
    np.random.seed(42)
    random_idx = np.sort(np.random.choice(adata_wt_full.n_obs, MAX_CELLS, replace=False))
    adata_wt = adata_wt_full[random_idx].to_memory()
else:
    print("Loading full WT dataset into memory...")
    adata_wt = adata_wt_full.to_memory()

try: adata_wt_full.file.close()
except: pass

# ==========================================
# 3. Preprocessing & Reference Mapping (Ingest)
# ==========================================
print("\n--- 3. Preprocessing & Reference Mapping ---")

# 1. Normalize and Log1p the Enriched data independently so it matches WT scale
print("Normalizing Enriched dataset...")
if adata_enr.X.max() > 100: 
    sc.pp.normalize_total(adata_enr, target_sum=1e4)
    sc.pp.log1p(adata_enr)

# 2. Align Features for Projection (Zero-Padding to bypass scanpy bugs)
print("Aligning features for Reference Mapping (Zero-padding missing genes)...")
adata_wt_empty = adata_wt[:0, :].copy()
adata_enr_aligned = ad.concat([adata_wt_empty, adata_enr], join='outer')
adata_enr_aligned = adata_enr_aligned[:, adata_wt.var_names].copy()

if not sp.issparse(adata_enr_aligned.X):
    adata_enr_aligned.X = sp.csr_matrix(adata_enr_aligned.X)

# 3. THE FIX: Forcefully Rebuild Missing Scanpy Metadata
print("Patching Reference Metadata for Ingest...")
if 'pca' not in adata_wt.uns:
    adata_wt.uns['pca'] = {}
adata_wt.uns['pca']['use_highly_variable'] = 'highly_variable' in adata_wt.var.columns

# If the PCA matrix itself was dropped from the file to save space, Ingest will fail.
# We recalculate the PCA/Neighbors here so Ingest can map the cells, 
# but we STRICTLY DO NOT run sc.tl.umap() so your original coordinates are perfectly preserved!
if 'PCs' not in adata_wt.varm or 'X_pca' not in adata_wt.obsm:
    print(" -> PCA loadings missing from WT file. Recomputing PCA & Neighbors safely...")
    sc.tl.pca(adata_wt, svd_solver='arpack', use_highly_variable='highly_variable' in adata_wt.var.columns)
    sc.pp.neighbors(adata_wt)

# 4. Project Enriched data directly onto the pristine WT Reference
print(" -> Running sc.tl.ingest (Projecting Enriched cells onto WT Manifold)...")
sc.tl.ingest(adata_enr_aligned, adata_wt, obs=CELL_TYPE_COL)

# 5. Extract the safely projected UMAP coordinates AND harmonized labels back to the main Enriched object
adata_enr.obsm['X_umap'] = adata_enr_aligned.obsm['X_umap']
adata_enr.obs[CELL_TYPE_COL] = adata_enr_aligned.obs[CELL_TYPE_COL]

del adata_wt_empty, adata_enr_aligned
import gc
gc.collect()

# ==========================================
# 3b. Concatenating Anchored Datasets
# ==========================================
print("\n--- 3b. Concatenating Datasets ---")
adata_wt.obs['dataset'] = 'WT'
adata_enr.obs['dataset'] = 'Enriched'

adata_wt.obs_names = adata_wt.obs_names + '_WT'
adata_enr.obs_names = adata_enr.obs_names + '_ENR'

adata = ad.concat([adata_wt, adata_enr], label='dataset_concat', keys=['WT', 'Enriched'], join='outer')

# Force the pristine, perfectly projected UMAP coordinates into the final object
adata.obsm['X_umap'] = np.vstack([adata_wt.obsm['X_umap'], adata_enr.obsm['X_umap']])
adata.raw = adata

del adata_wt, adata_enr
gc.collect()

print(f"Total cells after concatenation: {adata.n_obs:,}")

print("Saving new anchored UMAP plots...")
sc.pl.umap(adata, color='dataset', title='Dataset Origin (WT Reference vs Enriched Projected)', show=False)
plt.savefig("3b_UMAP_Dataset_Comparison.png", dpi=300, bbox_inches='tight')
plt.close()

# ==========================================
# 4. Global Subpopulation Differentials
# ==========================================
print(f"\n--- 4. Differential Expression across target types ---")
adata.obs['celltype_dataset'] = adata.obs[CELL_TYPE_COL].astype(str) + "_" + adata.obs['dataset'].astype(str)
adata.obs['celltype_dataset'] = adata.obs['celltype_dataset'].astype('category')

for cell_type in TARGET_CELL_TYPES:
    print(f"\nTesting: {cell_type}")
    group_enriched = f"{cell_type}_Enriched"
    group_wt = f"{cell_type}_WT"

    if group_enriched in adata.obs['celltype_dataset'].values and group_wt in adata.obs['celltype_dataset'].values:
        try:
            mask = adata.obs['celltype_dataset'].isin([group_enriched, group_wt])
            adata_sub = adata[mask].copy()
            
            # Downsample to prevent memory crash
            if adata_sub.n_obs > 10000:
                sc.pp.subsample(adata_sub, n_obs=10000, random_state=42)
                
            # Run optimized T-test
            sc.tl.rank_genes_groups(
                adata_sub, groupby='celltype_dataset', groups=[group_enriched], 
                reference=group_wt, method='t-test', use_raw=True, key_added=f'DE_{cell_type}'
            )
            
            adata.uns[f'DE_{cell_type}'] = adata_sub.uns[f'DE_{cell_type}']
            result = adata.uns[f'DE_{cell_type}']
            genes = result['names'][group_enriched]
            
            de_df = pd.DataFrame({
                'Gene': genes, 'Log2FC': result['logfoldchanges'][group_enriched],
                'P-val': result['pvals'][group_enriched], 'Adj_P-val': result['pvals_adj'][group_enriched]
            })
            
            # Use full dataset matrix for true accurate means
            mask_wt = adata.obs['celltype_dataset'] == group_wt
            mask_enr = adata.obs['celltype_dataset'] == group_enriched
            
            def get_mean_expr(matrix): return np.asarray(matrix.mean(axis=0)).flatten()
            
            de_df['Mean_Expr_WT'] = get_mean_expr(adata.raw[mask_wt, genes].X)
            de_df['Mean_Expr_Enriched'] = get_mean_expr(adata.raw[mask_enr, genes].X)
            
            de_df = de_df[['Gene', 'Mean_Expr_WT', 'Mean_Expr_Enriched', 'Log2FC', 'P-val', 'Adj_P-val']]
            safe_name = cell_type.replace(' ', '_').replace('+', 'pos').replace('/', '_')
            de_df.to_csv(f"4_DE_results_{safe_name}.csv", index=False)
            
            # Volcano Plot
            plt.figure(figsize=(8, 7))
            de_df['minuslog10pval'] = -np.log10(de_df['Adj_P-val'] + 1e-300)
            sig_up = (de_df['Adj_P-val'] < 0.05) & (de_df['Log2FC'] > 1)
            sig_dn = (de_df['Adj_P-val'] < 0.05) & (de_df['Log2FC'] < -1)
            
            plt.scatter(de_df.loc[~(sig_up | sig_dn), 'Log2FC'], de_df.loc[~(sig_up | sig_dn), 'minuslog10pval'], color='grey', alpha=0.3, s=10)
            plt.scatter(de_df.loc[sig_up, 'Log2FC'], de_df.loc[sig_up, 'minuslog10pval'], color='#d62728', alpha=0.7, s=15, label='Upregulated (Enriched)')
            plt.scatter(de_df.loc[sig_dn, 'Log2FC'], de_df.loc[sig_dn, 'minuslog10pval'], color='#1f77b4', alpha=0.7, s=15, label='Downregulated (WT)')
            
            viral_genes_df = de_df[de_df['Gene'].str.split('_').str[0].isin(KNOWN_EBV_GENES) | de_df['Gene'].str.startswith('EBV-') | de_df['Gene'].str.contains('TYPE1', case=False)]
            if not viral_genes_df.empty:
                plt.scatter(viral_genes_df['Log2FC'], viral_genes_df['minuslog10pval'], color='orange', edgecolor='black', s=40, label='EBV Genes', zorder=5)
                for _, row in viral_genes_df.sort_values('minuslog10pval', ascending=False).head(5).iterrows():
                    plt.text(row['Log2FC'] + 0.1, row['minuslog10pval'], row['Gene'], fontsize=9, fontweight='bold')
            
            plt.axvline(1, color='k', linestyle='--', linewidth=0.5)
            plt.axvline(-1, color='k', linestyle='--', linewidth=0.5)
            plt.axhline(-np.log10(0.05), color='k', linestyle='--', linewidth=0.5)
            plt.title(f"Volcano Plot: {cell_type} (Enriched vs WT)")
            plt.xlabel("Log2 Fold Change")
            plt.ylabel("-Log10(Adj P-val)")
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
            plt.tight_layout()
            plt.savefig(f"4_Volcano_Plot_{safe_name}.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            del adata_sub
            gc.collect()
            
        except Exception as e: print(f" -> Failed DE for {cell_type}. Error: {e}")

# ==========================================
# 5. Multi-Compartment Statistical Module
# ==========================================
print("\n--- 5. Running Multi-Compartment Statistical Module ---")

adata.obs['Disease_Group'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS' if 'Active' in str(x) or 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown')
)
adata.obs['Infection_Status'] = adata.obs.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV')
adata.obs['Day_Clean'] = adata.obs.get('Day', 'Unknown').astype(str).fillna('Unknown')

# Score ALL defined modules globally
for mod_name, genes in MODULE_GENES.items():
    valid_genes = [g for g in genes if g in adata.var_names]
    if valid_genes:
        sc.tl.score_genes(adata, valid_genes, score_name=mod_name, use_raw=False)
    else:
        adata.obs[mod_name] = 0.0

def analyze_compartment(adata, compartment_classes, compartment_name, primary_score):
    print(f"\n>> Analyzing Compartment: {compartment_name} (n={len(compartment_classes)} classes)")
    mask = (adata.obs['dataset'] == 'Enriched') & (adata.obs[CELL_TYPE_COL].isin(compartment_classes))
    adata_sub = adata[mask].copy()
    
    if adata_sub.n_obs < 10:
        print(f"Skipping {compartment_name}: Not enough cells ({adata_sub.n_obs}).")
        return

    ebv_genes_sub = [g for g in adata_sub.var_names if g.split('_')[0] in KNOWN_EBV_GENES or 'TYPE1' in g.upper() or 'EBV' in g.upper()]
    adata_sub.obs['total_ebv_counts'] = adata_sub[:, ebv_genes_sub].X.sum(axis=1) if ebv_genes_sub else 0

    print(f"Calculating Spearman correlations for {compartment_name}...")
    corr_results = []
    
    all_target_genes = []
    for g_list in MODULE_GENES.values(): all_target_genes.extend(g_list)
    valid_targets = [g for g in set(all_target_genes) if g in adata_sub.var_names]

    for disease in ['MS', 'HC']:
        df_disease = adata_sub[adata_sub.obs['Disease_Group'] == disease]
        viral_loads = df_disease.obs['total_ebv_counts'].values
        
        if df_disease.n_obs > 10 and np.std(viral_loads) > 0: 
            for gene in valid_targets:
                gene_expr = df_disease[:, gene].X.toarray().flatten() if sp.issparse(df_disease.X) else np.asarray(df_disease[:, gene].X).flatten()
                if np.std(gene_expr) > 0: 
                    rho, pval = stats.spearmanr(gene_expr, viral_loads)
                    corr_results.append({'Compartment': compartment_name, 'Disease_Group': disease, 'Gene': gene, 'Spearman_rho': rho, 'P-value': pval})

    if corr_results: 
        pd.DataFrame(corr_results).to_csv(f"5_Spearman_Corr_{compartment_name}_vs_ViralLoad.csv", index=False)
        print(f"Saved {compartment_name} correlation table.")

analyze_compartment(adata, B_CELL_CLASSES, "B_Cells", "ABC_Score")
analyze_compartment(adata, T_CELL_CLASSES, "T_Cells", "T_Exhaustion_Score")
analyze_compartment(adata, MONO_MAC_CLASSES, "Monocytes_Macrophages", "Mono_Inflammatory_Score")

# ==========================================
# 6. Save Final Anchored Object
# ==========================================
print("\n--- 6. Saving Final AnnData Object ---")

# ---> NEW: Drop any lingering unclassifiable cells before saving <---
adata = adata[adata.obs[CELL_TYPE_COL] != 'Unknown'].copy()

SAVE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
adata.write_h5ad(SAVE_PATH)
print(f"Pipeline complete. Full combined dataset saved to: {SAVE_PATH}")