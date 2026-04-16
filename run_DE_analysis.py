"""
Differential Expression Analysis Pipeline with Functional Enrichment

Performs pseudobulk-level differential expression analysis across multiple cell types
and biological contrasts (standard comparisons and statistical interactions). Includes
automated functional pathway enrichment for significant genes.

Outputs:
- Volcano plots (HTML + interactive)
- Significant genes (CSV)
- Enrichment results (HTML + CSV)
"""

import scanpy as sc
import decoupler as dc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import os
import plotly.express as px
import plotly.graph_objects as go
import gseapy as gp
import time
from datetime import datetime
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import warnings
warnings.filterwarnings('ignore')

# ==========================================
# 0. CONFIGURATION
# ==========================================
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
DIR_DE = "Dashboards_07_Differential_Expression_Pseudobulk"
DIR_DE_ENRICHMENT = "Dashboards_07_Differential_Expression_Pseudobulk/Functional_Enrichment"
os.makedirs(DIR_DE, exist_ok=True)
os.makedirs(DIR_DE_ENRICHMENT, exist_ok=True)

BIOLOGICAL_REPLICATE_COL = "PatientID"

# Validation thresholds
MIN_CELLS_PER_GROUP = 50
MIN_DONORS_PER_GROUP = 2
DE_PADJ_THRESHOLD = 0.05
DE_LOGFC_THRESHOLD = 1.0 

# Target Cell Types
TARGET_CELLS = {
    'CD16pos_NK': 'CD16+ NK cells',
    'ABC': 'Age-associated B cells',
    'NonClassMono': 'Non-classical monocytes',
    'Naive_Helper_T': 'Tcm/Naive helper T cells'
}

# --- LIST 1: STANDARD A vs B CONTRASTS ---
# Format: (Test_Group, Ref_Group, File_Prefix, Plot_Title)
STANDARD_CONTRASTS = [
    # Baseline Disease Signatures (Day 1, Mock)
    ('MS Active_Mock_Day 1', 'HC_Mock_Day 1', 'Base_Active_v_HC', 'Active MS vs HC (Day 1 Mock)'),
    ('MS Stable_Mock_Day 1', 'HC_Mock_Day 1', 'Base_Stable_v_HC', 'Stable MS vs HC (Day 1 Mock)'),
    
    # Direct Viral Reaction (EBV vs Mock per Day/Disease)
    ('MS Active_EBV_Day 7', 'MS Active_Mock_Day 7', 'Viral_Active_D7', 'Active MS: EBV vs Mock (Day 7)'),
    ('HC_EBV_Day 7', 'HC_Mock_Day 7', 'Viral_HC_D7', 'HC: EBV vs Mock (Day 7)'),
    ('MS Active_EBV_Day 15', 'MS Active_Mock_Day 15', 'Viral_Active_D15', 'Active MS: EBV vs Mock (Day 15)'),
    ('MS Stable_EBV_Day 15', 'MS Stable_Mock_Day 15', 'Viral_Stable_D15', 'Stable MS: EBV vs Mock (Day 15)'),
    ('HC_EBV_Day 15', 'HC_Mock_Day 15', 'Viral_HC_D15', 'HC: EBV vs Mock (Day 15)'),

    # Viral Kinetics (Tracking progression over time in infected cells)
    ('MS Active_EBV_Day 7', 'MS Active_EBV_Day 1', 'Kin_Active_D7_v_D1', 'Active MS: EBV Day 7 vs Day 1'),
    ('MS Active_EBV_Day 15', 'MS Active_EBV_Day 7', 'Kin_Active_D15_v_D7', 'Active MS: EBV Day 15 vs Day 7'),
    ('HC_EBV_Day 15', 'HC_EBV_Day 7', 'Kin_HC_D15_v_D7', 'HC: EBV Day 15 vs Day 7')
]

# --- LIST 2: STATISTICAL INTERACTION CONTRASTS ---
# Format: (Test_Disease, Ref_Disease, Target_Day, File_Prefix, Plot_Title)
# Question: "Does the Test Disease react differently to EBV than the Ref Disease at this Day?"
INTERACTION_CONTRASTS = [
    ('MS Active', 'HC', 'Day 7', 'Int_Active_v_HC_D7', 'Interaction: Active MS vs HC Viral Response (Day 7)'),
    ('MS Active', 'HC', 'Day 15', 'Int_Active_v_HC_D15', 'Interaction: Active MS vs HC Viral Response (Day 15)'),
    ('MS Stable', 'HC', 'Day 15', 'Int_Stable_v_HC_D15', 'Interaction: Stable MS vs HC Viral Response (Day 15)'),
    ('MS Active', 'MS Stable', 'Day 15', 'Int_Active_v_Stable_D15', 'Interaction: Active vs Stable MS Viral Response (Day 15)')
]

# Enrichment Analysis Settings
GENE_SETS = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
ENRICHMENT_PADJ_THRESHOLD = 0.05
ENRICHMENT_LOGFC_THRESHOLD = 0.5
ENRICHMENT_MIN_GENES = 10
ENRICHMENT_MAX_RETRIES = 2
ENRICHMENT_RETRY_WAIT = 2

# Logging/tracking
START_TIME = datetime.now()
RUN_STATS = {'contrasts_run': 0, 'contrasts_skipped': 0, 'enrichments_completed': 0}


# ==========================================
# 1. LOGGING HELPER
# ==========================================
def log_message(level, msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {msg}")

# Logging/tracking
START_TIME = datetime.now()
RUN_STATS = {'contrasts_run': 0, 'contrasts_skipped': 0, 'enrichments_completed': 0}


# ==========================================
# 1. LOGGING HELPER
# ==========================================
def log_message(level, msg):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{level}] {msg}")


# ==========================================
# 2. DATA VALIDATION HELPERS
# ==========================================
def validate_input_file(filepath):
    """Validate that input h5ad file exists and is readable."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
    try:
        adata_test = sc.read_h5ad(filepath)
        return True
    except Exception as e:
        raise ValueError(f"Cannot read h5ad file: {e}")


def validate_metadata_columns(adata, required_cols):
    """Check that required metadata columns exist."""
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required metadata columns: {missing}")
    return True


# ==========================================
# 3. DATA PREPARATION HELPER
# ==========================================
def prepare_counts_for_de(adata_ct):
    """Extract and validate count matrix for DE analysis."""
    if adata_ct.raw is not None:
        adata_ct = adata_ct.raw.to_adata()
    elif 'counts' in adata_ct.layers:
        adata_ct.X = adata_ct.layers['counts'].copy()
    
    # Ensure integer counts
    if sp.issparse(adata_ct.X):
        adata_ct.X.data = np.round(adata_ct.X.data).astype(int)
    else:
        adata_ct.X = np.round(adata_ct.X).astype(int)
    
    return adata_ct


# ==========================================
# 4. HTML PLOTTING HELPER
# ==========================================
def export_plot_and_table(fig, table_df, directory, filename, title):
    fig.update_layout(height=700, autosize=True, plot_bgcolor='white', margin=dict(l=20, r=20, t=50, b=50))
    fig_html = fig.to_html(full_html=False, include_plotlyjs='cdn', default_width='100%')
    
    table_html = f"""
    <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 40px; border-bottom: 2px solid #dee2e6; padding-bottom: 10px;">
        <h3 style="margin: 0;">Significant Results (padj < 0.05)</h3>
    </div>
    <div style="max-height: 500px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 4px;">
        {table_df.to_html(classes='data-table', index=False, float_format='%.4f')}
    </div>
    <style>
        body {{ font-family: -apple-system, sans-serif; padding: 20px; background: #f8f9fa; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); }}
        .data-table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
        .data-table th, .data-table td {{ padding: 12px; border-bottom: 1px solid #dee2e6; text-align: left; }}
        .data-table th {{ background-color: #e9ecef; position: sticky; top: 0; }}
    </style>
    """
    html_template = f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>{title}</title></head><body><div class='container'><h2>{title}</h2><div>{fig_html}</div>{table_html}</div></body></html>"
    
    with open(os.path.join(directory, filename), 'w', encoding='utf-8') as f:
        f.write(html_template)


# ==========================================
# 3. DATA VALIDATION HELPERS
# ==========================================
def validate_input_file(filepath):
    """Validate that input h5ad file exists and is readable."""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")
    try:
        adata_test = sc.read_h5ad(filepath)
        return True
    except Exception as e:
        raise ValueError(f"Cannot read h5ad file: {e}")


def validate_metadata_columns(adata, required_cols):
    """Check that required metadata columns exist."""
    missing = [c for c in required_cols if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required metadata columns: {missing}")
    return True


# ==========================================
# 4. DATA PREPARATION HELPER
# ==========================================
def prepare_counts_for_de(adata_ct):
    """Extract and validate count matrix for DE analysis."""
    if adata_ct.raw is not None:
        adata_ct = adata_ct.raw.to_adata()
    elif 'counts' in adata_ct.layers:
        adata_ct.X = adata_ct.layers['counts'].copy()
    
    # Ensure integer counts
    if sp.issparse(adata_ct.X):
        adata_ct.X.data = np.round(adata_ct.X.data).astype(int)
    else:
        adata_ct.X = np.round(adata_ct.X).astype(int)
    
    return adata_ct


# ==========================================
# 5. FUNCTIONAL ENRICHMENT HELPER
# ==========================================
def run_functional_enrichment(results_df, comparison_name, short_name, direction='upregulated'):
    """
    Run functional enrichment on significant genes from DE results.
    
    Parameters:
    - results_df: DE results dataframe with 'Gene', 'padj', 'log2FoldChange' columns
    - comparison_name: Name of the comparison (for file naming)
    - short_name: Cell type short name (CD16pos_NK, ABC, etc.)
    - direction: 'upregulated' or 'downregulated'
    
    Returns: enrichment results dataframe or None if insufficient genes
    """
    try:
        # Filter for significant genes based on direction
        if direction == 'upregulated':
            sig_genes = results_df[(results_df['padj'] < ENRICHMENT_PADJ_THRESHOLD) & 
                                   (results_df['log2FoldChange'] > ENRICHMENT_LOGFC_THRESHOLD)]['Gene'].tolist()
        else:
            sig_genes = results_df[(results_df['padj'] < ENRICHMENT_PADJ_THRESHOLD) & 
                                   (results_df['log2FoldChange'] < -ENRICHMENT_LOGFC_THRESHOLD)]['Gene'].tolist()
        
        # Filter out viral and artifact genes
        sig_genes = [g for g in sig_genes if not g.startswith('EBV-') and 'TYPE1' not in g.upper()]
        
        if len(sig_genes) < ENRICHMENT_MIN_GENES:
            print(f"         ⏭️ Enrichment skipped: Only {len(sig_genes)} genes (need ≥{ENRICHMENT_MIN_GENES})")
            return None
        
        print(f"         Running enrichment on {len(sig_genes)} {direction} genes...")
        
        # Run enrichment with retry logic
        df_sig = None
        for attempt in range(ENRICHMENT_MAX_RETRIES):
            try:
                # Run enrichment (gseapy manages its own timeouts)
                print(f"            [Attempt {attempt + 1}/{ENRICHMENT_MAX_RETRIES}]", flush=True)
                enr = gp.enrichr(gene_list=sig_genes, gene_sets=GENE_SETS, outdir=None)
                df_enr = enr.res2d
                df_sig = df_enr[df_enr['Adjusted P-value'] < ENRICHMENT_PADJ_THRESHOLD].sort_values('Adjusted P-value')
                
                if df_sig.empty:
                    print(f"         ℹ️ No significant pathways (padj < {ENRICHMENT_PADJ_THRESHOLD})")
                    return None
                
                # Successfully completed
                break
                
            except Exception as e:
                print(f"            Attempt {attempt + 1} failed: {str(e)[:100]}")
                if attempt < ENRICHMENT_MAX_RETRIES - 1:
                    time.sleep(ENRICHMENT_RETRY_WAIT)  # Wait before retry
                else:
                    print(f"         ⚠️ Enrichment failed after {ENRICHMENT_MAX_RETRIES} attempts")
                    return None
        
        if df_sig is None or df_sig.empty:
            return None
        
        # Clean up pathway names
        df_sig['Term_Clean'] = df_sig['Term'].apply(lambda x: x.split(' (GO:')[0] if '(GO:' in str(x) else str(x))
        
        # Save to CSV
        csv_filename = f"{short_name}_{comparison_name}_{direction}_Enrichment.csv"
        df_sig.to_csv(os.path.join(DIR_DE_ENRICHMENT, csv_filename), index=False)
        print(f"         ✅ Saved {len(df_sig)} significant pathways to {csv_filename}")
        
        return df_sig
        
    except Exception as e:
        print(f"         ⚠️ Enrichment failed: {e}")
        return None


# ==========================================
# 6. ENRICHMENT EXPORT HELPER
# ==========================================
def export_enrichment_plot_and_table(enr_df, comparison_name, short_name, direction, title):
    """
    Export enrichment results as interactive plot and HTML table.
    """
    if enr_df is None or enr_df.empty:
        return
    
    try:
        # Create volcano-style scatter for enrichment
        enr_df['-log10(Adj_P-val)'] = -np.log10(enr_df['Adjusted P-value'] + 1e-300)
        
        # Parse overlap (e.g., "5/120" -> numerator)
        enr_df['Overlap_Count'] = enr_df['Overlap'].str.split('/').str[0].astype(int)
        enr_df['Pathway_Size'] = enr_df['Overlap'].str.split('/').str[1].astype(int)
        
        # Create scatter plot
        fig = px.scatter(
            enr_df.head(50), 
            x='Overlap_Count', 
            y='-log10(Adj_P-val)', 
            size='Pathway_Size',
            hover_data=['Term_Clean', 'Gene_set', 'Adjusted P-value'],
            title=f"{short_name}: Functional Enrichment ({direction})",
            labels={'Overlap_Count': 'Overlapping Genes', '-log10(Adj_P-val)': '-log10(Adj P-value)'}
        )
        fig.add_hline(y=-np.log10(ENRICHMENT_PADJ_THRESHOLD), line_dash="dash", line_color="gray", 
                      annotation_text=f"padj={ENRICHMENT_PADJ_THRESHOLD}")
        
        # Annotate top 10 pathways
        top_pathways = enr_df.nsmallest(10, 'Adjusted P-value')
        for _, row in top_pathways.iterrows():
            fig.add_annotation(
                x=row['Overlap_Count'], 
                y=row['-log10(Adj_P-val)'],
                text=row['Term_Clean'][:30] + ('...' if len(row['Term_Clean']) > 30 else ''),
                showarrow=False,
                yshift=10,
                font=dict(size=9)
            )
        
        # Export
        filename = f"{short_name}_{comparison_name}_{direction}_Enrichment.html"
        export_plot_and_table(fig, enr_df[['Term_Clean', 'Gene_set', 'Adjusted P-value', 'Overlap', 'Genes']].head(200), 
                             DIR_DE_ENRICHMENT, filename, f"{short_name}: {direction} Enrichment - {title}")
        
    except Exception as e:
        print(f"         ⚠️ Plot export failed: {e}")


# ==========================================
# 7. LOAD & PREP DATA
# ==========================================
log_message("INFO", f"Loading data from {FILE_PATH}...")
try:
    validate_input_file(FILE_PATH)
except (FileNotFoundError, ValueError) as e:
    log_message("ERROR", str(e))
    exit(1)

adata = sc.read_h5ad(FILE_PATH)
log_message("INFO", f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

# Filter to WT dataset
if 'dataset' not in adata.obs.columns:
    log_message("ERROR", "'dataset' column not found in metadata")
    exit(1)

adata = adata[adata.obs['dataset'] == 'WT'].copy()
log_message("INFO", f"Retained {adata.n_obs:,} WT cells")

adata.obs[BIOLOGICAL_REPLICATE_COL] = adata.obs[BIOLOGICAL_REPLICATE_COL].astype(str)

adata.obs['Disease_Group'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else 'HC')
).astype('object')
adata.obs['Infection_Clean'] = adata.obs.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV').astype('object')
adata.obs['Day_Clean'] = adata.obs['Day'].astype(str).astype('object')
adata.obs['DE_Combo'] = adata.obs['Disease_Group'] + "_" + adata.obs['Infection_Clean'] + "_" + adata.obs['Day_Clean']


# ==========================================
# 8. ENGINE A: STANDARD CONTRASTS
# ==========================================
def run_standard_contrast(adata_ct, test_group, ref_group, prefix, title, short_name):
    print(f"\n   [Standard] Testing: {test_group} vs {ref_group}")
    adata_test = adata_ct[adata_ct.obs['DE_Combo'].isin([test_group, ref_group])].copy()
    
    if adata_test.n_obs < MIN_CELLS_PER_GROUP:
        print(f"      ⏭️ SKIPPED: Insufficient cells ({adata_test.n_obs} < {MIN_CELLS_PER_GROUP}).")
        RUN_STATS['contrasts_skipped'] += 1
        return

    t_donors = adata_test.obs[adata_test.obs['DE_Combo'] == test_group][BIOLOGICAL_REPLICATE_COL].nunique()
    r_donors = adata_test.obs[adata_test.obs['DE_Combo'] == ref_group][BIOLOGICAL_REPLICATE_COL].nunique()
    if t_donors < MIN_DONORS_PER_GROUP or r_donors < MIN_DONORS_PER_GROUP:
        print(f"      ⏭️ SKIPPED: Insufficient donors (Test: {t_donors}, Ref: {r_donors}, need ≥{MIN_DONORS_PER_GROUP}per group).")
        RUN_STATS['contrasts_skipped'] += 1
        return

    try:
        pbulk = dc.pp.pseudobulk(adata_test, sample_col=BIOLOGICAL_REPLICATE_COL, groups_col='DE_Combo', mode='sum')
        counts = pd.DataFrame(pbulk.X.toarray() if sp.issparse(pbulk.X) else pbulk.X, index=pbulk.obs_names, columns=pbulk.var_names)
        metadata = pbulk.obs[['DE_Combo', BIOLOGICAL_REPLICATE_COL]].copy()

        dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors='DE_Combo', quiet=True)
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=['DE_Combo', test_group, ref_group], quiet=True)
        stat_res.summary()
        process_and_export_results(stat_res.results_df, test_group, ref_group, prefix, title, short_name)
        RUN_STATS['contrasts_run'] += 1
    except Exception as e:
        log_message("ERROR", f"[Standard] {test_group} vs {ref_group}: {e}")
        RUN_STATS['contrasts_skipped'] += 1


# ==========================================
# 9. ENGINE B: STATISTICAL INTERACTION
# ==========================================
def run_interaction_contrast(adata_ct, test_disease, ref_disease, target_day, prefix, title, short_name):
    print(f"\n   [Interaction] Testing: {test_disease} vs {ref_disease} response at {target_day}")
    
    # Subset to the specific day and the two diseases (keeping both Mock and EBV)
    mask = (adata_ct.obs['Day_Clean'] == target_day) & (adata_ct.obs['Disease_Group'].isin([test_disease, ref_disease]))
    adata_test = adata_ct[mask].copy()

    if adata_test.n_obs < MIN_CELLS_PER_GROUP:
        print(f"      ⏭️ SKIPPED: Insufficient cells ({adata_test.n_obs} < {MIN_CELLS_PER_GROUP}).")
        RUN_STATS['contrasts_skipped'] += 1
        return

    # Check that we have enough donors in all 4 required quadrants
    quadrants = adata_test.obs.groupby(['Disease_Group', 'Infection_Clean'])[BIOLOGICAL_REPLICATE_COL].nunique()
    min_donors_ok = all(quadrants.get((d, i), 0) >= MIN_DONORS_PER_GROUP for d in [test_disease, ref_disease] for i in ['Mock', 'EBV'])
    if not min_donors_ok:
        print(f"      ⏭️ SKIPPED: Insufficient donors across interaction quadrants.")
        print(f"         Available donors: \n{quadrants}")
        RUN_STATS['contrasts_skipped'] += 1
        return

    try:
        # Create a combined grouping column just for pseudo-bulking aggregation
        adata_test.obs['PB_Group'] = adata_test.obs['Disease_Group'] + "_" + adata_test.obs['Infection_Clean']
        
        pbulk = dc.pp.pseudobulk(adata_test, sample_col=BIOLOGICAL_REPLICATE_COL, groups_col='PB_Group', mode='sum')
        counts = pd.DataFrame(pbulk.X.toarray() if sp.issparse(pbulk.X) else pbulk.X, index=pbulk.obs_names, columns=pbulk.var_names)
        
        # Build clean metadata for the interaction formula
        metadata = pbulk.obs[['Disease_Group', 'Infection_Clean']].copy()
        
        # Explicitly set the reference levels so the math subtracts in the correct direction
        metadata['Disease_Group'] = pd.Categorical(metadata['Disease_Group'], categories=[ref_disease, test_disease])
        metadata['Infection_Clean'] = pd.Categorical(metadata['Infection_Clean'], categories=['Mock', 'EBV'])

        # Run the interaction formula: ~ Disease + Infection + (Disease x Infection)
        dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="Disease_Group * Infection_Clean", quiet=True)
        dds.deseq2()

        # Dynamically find the interaction coefficient column generated by Patsy
        interaction_coef = [c for c in dds.obsm['design_matrix'].columns if ':' in c][0]
        
        # Extract stats using the specific interaction coefficient name instead of a contrast list
        stat_res = DeseqStats(dds, name=interaction_coef, quiet=True)
        stat_res.summary()
        
        # Define labels for the volcano plot
        up_label = f"Hyper-responsive in {test_disease}"
        down_label = f"Failed/Blunted in {test_disease}"
        process_and_export_results(stat_res.results_df, up_label, down_label, prefix, title, short_name)
        RUN_STATS['contrasts_run'] += 1

    except Exception as e:
        log_message("ERROR", f"[Interaction] {test_disease} vs {ref_disease} at {target_day}: {e}")
        RUN_STATS['contrasts_skipped'] += 1


# ==========================================
# 10. RESULT FORMATTING & EXPORT
# ==========================================
def process_and_export_results(results_df, up_label, down_label, prefix, title, short_name):
    res_df = results_df.dropna(subset=['padj']).reset_index()
    res_df.rename(columns={'index': 'Gene'}, inplace=True)
    
    res_df['-log10(Adj_P-val)'] = -np.log10(res_df['padj'] + 1e-300)
    res_df['Significance'] = 'Not Sig'
    
    # Apply significance filter using configured thresholds
    sig_mask = (res_df['padj'] < DE_PADJ_THRESHOLD) & (abs(res_df['log2FoldChange']) > DE_LOGFC_THRESHOLD)
    res_df.loc[sig_mask & (res_df['log2FoldChange'] > 0), 'Significance'] = up_label
    res_df.loc[sig_mask & (res_df['log2FoldChange'] < 0), 'Significance'] = down_label
    
    sig_genes_df = res_df[res_df['Significance'] != 'Not Sig'].sort_values('padj')
    n_sig = len(sig_genes_df)
    
    # Create volcano plot
    fig_volc = px.scatter(
        res_df, x='log2FoldChange', y='-log10(Adj_P-val)', color='Significance', 
        hover_data=['Gene', 'padj', 'baseMean'],
        color_discrete_map={'Not Sig': '#e0e0e0', up_label: '#d62728', down_label: '#1f77b4'},
        title=f"{short_name}: {title}"
    )
    
    # Annotate top genes
    top_genes = sig_genes_df.head(15)
    for _, row in top_genes.iterrows():
        fig_volc.add_annotation(x=row['log2FoldChange'], y=row['-log10(Adj_P-val)'], text=row['Gene'], showarrow=False, yshift=10)
        
    filename = f"{short_name}_{prefix}_Volcano.html"
    export_plot_and_table(fig_volc, sig_genes_df.head(200), DIR_DE, filename, f"{short_name} DE: {title}")
    
    csv_name = f"{short_name}_{prefix}_Sig_Genes.csv"
    sig_genes_df.to_csv(os.path.join(DIR_DE, csv_name), index=False)
    print(f"      ✅ Success. Found {n_sig} sig genes ({len(top_genes)} annotated) in {filename}")
    
    # ==========================================
    # FUNCTIONAL ENRICHMENT ANALYSIS
    # ==========================================
    try:
        log_message("INFO", f"Running enrichment for {short_name} - {prefix}")
        
        # Upregulated genes enrichment
        enr_up = run_functional_enrichment(res_df, prefix, short_name, direction='upregulated')
        if enr_up is not None:
            export_enrichment_plot_and_table(enr_up, prefix, short_name, 'upregulated', title)
            RUN_STATS['enrichments_completed'] += 1
            time.sleep(1.0)  # Rate limiting for gseapy API
        
        # Downregulated genes enrichment
        enr_down = run_functional_enrichment(res_df, prefix, short_name, direction='downregulated')
        if enr_down is not None:
            export_enrichment_plot_and_table(enr_down, prefix, short_name, 'downregulated', title)
            RUN_STATS['enrichments_completed'] += 1
            time.sleep(1.0)  # Rate limiting for gseapy API
            
    except Exception as e:
        log_message("WARNING", f"Enrichment failed (continuing with DE): {e}")


# ==========================================
# 11. EXECUTION LOOP
# ==========================================
log_message("INFO", "="*50)
log_message("INFO", "🚀 INITIATING MASTER DE PIPELINE")
log_message("INFO", "="*50)

for short_name, exact_cell_type in TARGET_CELLS.items():
    log_message("INFO", f"\n[{short_name.upper()}] Isolating {exact_cell_type}...")
    
    # Filter to target cell type
    if 'majority_voting' not in adata.obs.columns:
        log_message("ERROR", "'majority_voting' column not found")
        continue
    
    adata_ct = adata[adata.obs['majority_voting'] == exact_cell_type].copy()
    if adata_ct.n_obs == 0:
        log_message("WARNING", f"No cells found for {exact_cell_type}")
        continue
    
    log_message("INFO", f"Found {adata_ct.n_obs:,} cells for {exact_cell_type}")
    
    # Prepare counts
    adata_ct = prepare_counts_for_de(adata_ct)
    
    # Run Standard Comparisons
    for test_grp, ref_grp, prefix, title in STANDARD_CONTRASTS:
        run_standard_contrast(adata_ct, test_grp, ref_grp, prefix, title, short_name)

    # Run Statistical Interactions
    for test_dis, ref_dis, target_day, prefix, title in INTERACTION_CONTRASTS:
        run_interaction_contrast(adata_ct, test_dis, ref_dis, target_day, prefix, title, short_name)

# Print summary statistics
elapsed_time = (datetime.now() - START_TIME).total_seconds() / 60
log_message("INFO", "="*50)
log_message("INFO", "🎉 PIPELINE COMPLETE")
log_message("INFO", f"Results saved to: {DIR_DE}/")
log_message("INFO", f"\n📊 SUMMARY STATISTICS:")
log_message("INFO", f"   Contrasts Run: {RUN_STATS['contrasts_run']}")
log_message("INFO", f"   Contrasts Skipped: {RUN_STATS['contrasts_skipped']}")
log_message("INFO", f"   Enrichments Completed: {RUN_STATS['enrichments_completed']}")
log_message("INFO", f"   Total Runtime: {elapsed_time:.1f} minutes")
log_message("INFO", "="*50)