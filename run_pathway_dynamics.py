import scanpy as sc
import pandas as pd
import numpy as np
import gseapy as gp
import os
import time
import re
import warnings

warnings.filterwarnings('ignore')

# ==========================================
# 1. Setup & Configuration
# ==========================================
print("Loading combined integrated h5ad file...")
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
adata = sc.read_h5ad(FILE_PATH)

# Clean var names just in case they haven't been saved yet
adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
adata.var_names_make_unique()

OUTPUT_DIR = "Dashboards_05_Pathway_Dynamics"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define Compartments
COMPARTMENTS = {
    'B-Cells': 'b cell|b|plasma|plasmablast',
    'T-Cells': 't cell|tcm|tem|treg|regulatory|cycling t',
    'Myeloid': 'monocyte|macrophage'
}

# Standardize Metadata
adata.obs['Disease_Group'] = adata.obs.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
)
adata.obs['Infection_Clean'] = adata.obs.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV')

def get_day_sort_key(day_str):
    match = re.search(r'\d+', str(day_str))
    return int(match.group()) if match else 9999

adata.obs['Day_Sort'] = adata.obs.get('Day', 'Unknown').astype(str).apply(get_day_sort_key)

TARGET_DAYS_MAP = {1: 'Day 1', 7: 'Day 7', 15: 'Day 15'}
DISEASES = ['HC', 'MS Active', 'MS Stable']
INFECTIONS = ['Mock', 'EBV', 'Combined']

GENE_SETS = ['GO_Biological_Process_2023', 'KEGG_2021_Human']

# Isolate ONLY the Whole Transcriptome dataset
adata_wt = adata[adata.obs['dataset'] == 'WT'].copy()
print(f"Isolated WT dataset: {adata_wt.n_obs} total cells.")

# ==========================================
# 2. Extract Signatures & Query Enrichr
# ==========================================
master_results = []
html_summary = []

for comp_name, comp_regex in COMPARTMENTS.items():
    print(f"\n================ Processing {comp_name} ================")
    
    # Isolate the compartment within WT
    mask_comp = adata_wt.obs['majority_voting'].astype(str).str.contains(comp_regex, case=False, na=False)
    adata_comp = adata_wt[mask_comp].copy()
    
    if adata_comp.n_obs == 0:
        continue
        
    for day_int, day_str in TARGET_DAYS_MAP.items():
        for disease in DISEASES:
            for inf in INFECTIONS:
                
                strat_name = f"{day_str}-{disease}-{inf}"
                
                # Build the specific group mask
                m_day = adata_comp.obs['Day_Sort'] == day_int
                m_dis = adata_comp.obs['Disease_Group'] == disease
                
                if inf == 'Combined':
                    m_inf = adata_comp.obs['Infection_Clean'].isin(['Mock', 'EBV'])
                else:
                    m_inf = adata_comp.obs['Infection_Clean'] == inf
                    
                m_target = m_day & m_dis & m_inf
                n_cells = m_target.sum()
                
                # Require at least 20 cells to calculate a meaningful pathway signature
                if n_cells < 20:
                    print(f" -> Skipping {strat_name} (Only {n_cells} cells)")
                    continue
                    
                print(f" -> Analyzing: {strat_name} ({n_cells} cells)")
                
                # Create binary comparison: This Specific Stratification vs ALL OTHER cells in the compartment
                adata_comp.obs['Comparison'] = np.where(m_target, 'Target', 'Background')
                
                try:
                    sc.tl.rank_genes_groups(adata_comp, groupby='Comparison', groups=['Target'], reference='Background', method='t-test', use_raw=True)
                    
                    res = adata_comp.uns['rank_genes_groups']
                    genes = res['names']['Target']
                    pvals = res['pvals_adj']['Target']
                    logfcs = res['logfoldchanges']['Target']
                    
                    # Positive markers only
                    sig_mask = (pvals < 0.05) & (logfcs > 0.5)
                    top_genes = list(genes[sig_mask][:150])
                    
                    # Remove viral genes from host pathway queries
                    top_genes = [g for g in top_genes if not g.startswith('EBV-') and 'TYPE1' not in g.upper()]
                    num_input_genes = len(top_genes)
                    
                    if num_input_genes >= 10:
                        enr = gp.enrichr(gene_list=top_genes, gene_sets=GENE_SETS, outdir=None)
                        df_enr = enr.res2d
                        df_sig = df_enr[df_enr['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
                        
                        if not df_sig.empty:
                            df_sig['Term_Clean'] = df_sig['Term'].apply(lambda x: x.split(' (GO:')[0])
                            
                            pathway_html_list = []
                            genes_html_list = []
                            
                            for idx, row in df_sig.head(10).iterrows():
                                # Parse the "15/200" overlap string
                                overlap_str = str(row['Overlap'])
                                try:
                                    over_genes, total_path_genes = overlap_str.split('/')
                                except:
                                    over_genes, total_path_genes = "?", "?"
                                    
                                master_results.append({
                                    'Compartment': comp_name, 'Stratification': strat_name, 
                                    'Day': day_str, 'Disease': disease, 'Infection': inf,
                                    'Cells_in_Group': n_cells, 'Input_Genes': num_input_genes,
                                    'Database': row['Gene_set'], 'Pathway': row['Term_Clean'], 
                                    'Adj_P_Value': row['Adjusted P-value'], 
                                    'Pathway_Total_Genes': total_path_genes,
                                    'Overlapping_Genes': over_genes,
                                    'Genes_Matched': row['Genes']
                                })
                                
                                if idx < 3: # Keep top 3 for the dashboard
                                    # Format Pathway String
                                    pathway_html_list.append(
                                        f"• <b>{row['Term_Clean']}</b><br>"
                                        f"<i>&nbsp;&nbsp;p={row['Adjusted P-value']:.2e} | Overlap: {over_genes}/{total_path_genes}</i>"
                                    )
                                    
                                    # Extract and Format the Top 20 Genes string
                                    matched_genes_list = str(row['Genes']).split(';')
                                    top_20 = ", ".join(matched_genes_list[:20])
                                    if len(matched_genes_list) > 20:
                                        top_20 += "..."
                                        
                                    genes_html_list.append(
                                        f"<div style='margin-bottom: 12px;'>"
                                        f"<b>{row['Term_Clean']}</b>:<br>"
                                        f"<span style='font-size: 0.85em; color: #555; line-height: 1.2;'>{top_20}</span>"
                                        f"</div>"
                                    )
                                    
                            html_summary.append({
                                'Compartment': comp_name,
                                'Stratification': f"<b>{strat_name}</b>",
                                'Group Cells': f"{n_cells:,}",
                                'Input Marker Genes': num_input_genes,
                                'Top Active Pathways': "<br><br>".join(pathway_html_list),
                                'Top 20 Overlapping Genes': "".join(genes_html_list)
                            })
                        else:
                            html_summary.append({
                                'Compartment': comp_name, 'Stratification': f"<b>{strat_name}</b>",
                                'Group Cells': f"{n_cells:,}", 'Input Marker Genes': num_input_genes,
                                'Top Active Pathways': "<i>No significant pathways identified.</i>",
                                'Top 20 Overlapping Genes': "<i>N/A</i>"
                            })
                        time.sleep(1.0) # Rate limit
                    else:
                        html_summary.append({
                                'Compartment': comp_name, 'Stratification': f"<b>{strat_name}</b>",
                                'Group Cells': f"{n_cells:,}", 'Input Marker Genes': num_input_genes,
                                'Top Active Pathways': "<i>Not enough marker genes (<10) for pathway analysis.</i>",
                                'Top 20 Overlapping Genes': "<i>N/A</i>"
                            })
                        
                except Exception as e:
                    print(f"      Error testing {strat_name}: {e}")

# ==========================================
# 3. Export Results
# ==========================================
print("\n--- Saving Results ---")

df_master = pd.DataFrame(master_results)
if not df_master.empty:
    csv_path = os.path.join(OUTPUT_DIR, "Master_Pathway_Stratifications_WT.csv")
    df_master.to_csv(csv_path, index=False)
    print(f"Saved complete deep-stat data to {csv_path}")

df_html = pd.DataFrame(html_summary)
if not df_html.empty:
    html_table = df_html.to_html(escape=False, index=False, classes='styled-table')
    
    html_template = f"""
    <!DOCTYPE html>
    <html><head><meta charset="utf-8"><title>WT Stratified Pathway Dynamics</title>
    <style>
        body {{ font-family: -apple-system, sans-serif; margin: 0; padding: 30px; background-color: #f8f9fa; color: #333; }}
        h1 {{ text-align: center; color: #2c3e50; margin-bottom: 5px; }}
        p.subtitle {{ text-align: center; color: #6c757d; margin-bottom: 30px; }}
        .styled-table {{ border-collapse: collapse; margin: 25px 0; font-size: 0.9em; width: 100%; box-shadow: 0 0 20px rgba(0, 0, 0, 0.05); background-color: white; table-layout: fixed; }}
        .styled-table thead tr {{ background-color: #0d6efd; color: #ffffff; text-align: left; }}
        .styled-table th, .styled-table td {{ padding: 15px 20px; border-bottom: 1px solid #dddddd; vertical-align: top; line-height: 1.4; word-wrap: break-word; }}
        /* Set specific column widths to ensure the gene lists don't stretch the table too wide */
        .styled-table th:nth-child(5) {{ width: 30%; }} 
        .styled-table th:nth-child(6) {{ width: 35%; }}
        .styled-table tbody tr:nth-of-type(even) {{ background-color: #f8f9fa; }}
        .styled-table tbody tr:hover {{ background-color: #e7f1ff; }}
        td b {{ color: #212529; }}
        td i {{ color: #6c757d; font-size: 0.9em; }}
    </style></head>
    <body>
        <h1>🧬 WT Stratified Pathway Enrichment</h1>
        <p class="subtitle">Defining biological processes for specific stratifications vs background compartment</p>
        {html_table}
    </body></html>
    """
    
    html_path = os.path.join(OUTPUT_DIR, "Interactive_Pathway_Dynamics_WT.html")
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_template)
    print(f"Saved interactive dashboard to {html_path}")

print("Stratified Pathway analysis complete!")
