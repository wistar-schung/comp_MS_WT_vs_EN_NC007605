import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
import scipy.stats as stats
import scipy.sparse as sp
import gseapy as gp
import re
import os
import time
import warnings
import argparse

warnings.filterwarnings('ignore')

# ==========================================
# argparse setup
# ==========================================
parser = argparse.ArgumentParser(description="Generate Master Dashboard")
parser.add_argument("--section", type=str, default="all", 
                    help="Run only specific sections (e.g., '4', 'global', '4,5', 'all'). Sections: 4=Global, 5=Compartments, 6=EBV, 7=Day, 8=Targeted, 9=DE, 10=Pathway")
args = parser.parse_args()

sections_to_run = [s.strip().lower() for s in args.section.split(',')]

def should_run(sec_id, name=None):
    if "all" in sections_to_run:
        return True
    return str(sec_id) in sections_to_run or (name and name.lower() in sections_to_run)

# ==========================================
# 0. Setup Directories & Load Data
# ==========================================
print("Setting up output directories...")
DIR_GLOBAL = "Dashboards_01_Global"
DIR_COMPARTMENTS = "Dashboards_02_Compartments"
DIR_EBV = "Dashboards_03_EBV_DeepDive"
DIR_DAYS = "Dashboards_04_Day_Contrasts"
DIR_PATHWAYS = "Dashboards_05_Pathway_Dynamics"
DIR_TARGETS = "Dashboards_06_Targeted_Deep_Dives"
DIR_DE = "Dashboards_07_Differential_Expression"

for d in [DIR_GLOBAL, DIR_COMPARTMENTS, DIR_EBV, DIR_DAYS, DIR_PATHWAYS, DIR_TARGETS, DIR_DE]:
    os.makedirs(d, exist_ok=True)

print("Loading processed h5ad file...")
FILE_PATH = "/Users/schung/work/comp_MS_WT_vs_EN_NC007605/combined_integrated_analysis_qc_MULTI_compartment.h5ad"
adata = sc.read_h5ad(FILE_PATH)

# --- CLEAN UP VIRAL GENE NAMES GLOBALLY ---
print("Cleaning up _TYPE1 suffixes from gene names...")
adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
adata.var_names_make_unique()

# ==========================================
# 1. Targeted Gene Modules & Scoring
# ==========================================
ABC_GENES = ["ITGAX", "TBX21", "FCRL5", "ZEB2", "CXCR3", "LILRB1", "LILRB2"]
T_EXHAUSTION_GENES = ["PDCD1", "LAG3", "HAVCR2", "CTLA4", "TIGIT", "TOX", "GZMB", "PRF1", "IFNG"]
MONO_INFLAMMATORY_GENES = ["IL1B", "TNF", "CXCL8", "CCL2", "NLRP3", "TLR4", "S100A8", "S100A9"]
IFN_GENES = ["ISG15", "IFI44L", "IFIT1", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "IRF7", "STAT1", "BST2", "XAF1"] 

KNOWN_EBV_GENES = ["A73", "BALF1", "BALF2", "BALF3", "BALF4", "BALF5", "BARF1", "BBLF1", "BBLF2/BBLF3", "BBLF4", "BBRF1", "BBRF2", "BBRF3", "BCRF1", "BDLF1", "BDLF2", "BDLF3", "BDLF3.5", "BDLF4", "BFLF1", "BFLF2", "BFRF1", "BFRF1A", "BFRF2", "BFRF3", "BGLF1", "BGLF2", "BGLF3", "BGLF3.5", "BGLF4", "BGLF5", "BGRF1/BDRF1", "BHRF1", "BILF1", "BILF2", "BKRF2", "BKRF3", "BKRF4", "BLLF1", "BLLF2", "BLLF3", "BLRF1", "BLRF2", "BMRF1", "BMRF2", "BNLF2a", "BNLF2b", "BNRF1", "BORF1", "BORF2", "BPLF1", "BRLF1", "BRRF1", "BRRF2", "BSLF1", "BSLF2/BMLF1", "BSRF1", "BTRF1", "BVLF1", "BVRF1", "BVRF2", "BXLF1", "BXLF2", "BXRF1", "BZLF1", "BZLF2", "BaRF1", "BcLF1", "BcRF1", "BdRF1", "EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B/EBNA-3C", "EBNA-LP", "LF1", "LF2", "LMP-1", "LMP-2A", "LMP-2B", "RPMS1"]
LATENT_GENES = ['LMP-1', 'LMP-2A', 'LMP-2B', 'EBNA-1', 'EBNA-2', 'EBNA-3A', 'EBNA-3B', 'EBNA-3C', 'EBNA-LP']
LYTIC_GENES  = ['BZLF1', 'BZLF2', 'LF1', 'LF3']

print("Calculating compartment scores and pure viral load on the fly...")
def score_if_valid(adata, genes, name):
    valid = [g for g in genes if g in adata.var_names]
    if valid: sc.tl.score_genes(adata, valid, score_name=name, use_raw=False)
    else: adata.obs[name] = 0.0

score_if_valid(adata, ABC_GENES, 'ABC_Score')
score_if_valid(adata, T_EXHAUSTION_GENES, 'T_Exhaustion_Score')
score_if_valid(adata, MONO_INFLAMMATORY_GENES, 'Mono_Inflammatory_Score')
score_if_valid(adata, IFN_GENES, 'IFN_Score')
score_if_valid(adata, LATENT_GENES, 'Latent_Score')
score_if_valid(adata, LYTIC_GENES, 'Lytic_Score')

ebv_genes_in_data = [g for g in adata.var_names if g in KNOWN_EBV_GENES or g.startswith('EBV-')]

ebv_genes = [g for g in adata.var_names if g in KNOWN_EBV_GENES or g.startswith('EBV-')]
if ebv_genes:
    adata.obs['Plot_Viral_Load'] = np.ravel(adata[:, ebv_genes].X.sum(axis=1))
else:
    adata.obs['Plot_Viral_Load'] = 0

adata.obs['Plot_Viral_Load_Log'] = np.log1p(adata.obs['Plot_Viral_Load'])

def assign_state(row):
    if row.get('Plot_Viral_Load', 0) < 2: return "Negative"
    elif row.get('Lytic_Score', 0) > row.get('Latent_Score', 0): return "Lytic-Like"
    else: return "Latent-Like"
adata.obs['Viral_State'] = adata.obs.apply(assign_state, axis=1)

# ==========================================
# 2. Build Master Plotting DF + Strict Parse Intersection
# ==========================================
print("Extracting metadata for plotting...")
plot_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)

cols_to_copy = [
    'majority_voting', 'dataset', 'Disease_Condition (Detail)', 'Day', 'Infection', 
    'ABC_Score', 'T_Exhaustion_Score', 'Mono_Inflammatory_Score', 'IFN_Score', 
    'Plot_Viral_Load', 'Plot_Viral_Load_Log', 'Latent_Score', 'Lytic_Score', 'Viral_State',
    'conf_score'
]

for c in cols_to_copy:
    if c in adata.obs.columns: plot_df[c] = adata.obs[c].values

plot_df['Disease_Group'] = plot_df.get('Disease_Condition (Detail)', 'Unknown').apply(
    lambda x: 'MS Active' if 'Active' in str(x) else ('MS Stable' if 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else 'Unknown'))
)
plot_df['Infection_Clean'] = plot_df.get('Infection', 'Unknown').apply(lambda x: 'Mock' if 'Mock' in str(x) else 'EBV')
plot_df['Day_Clean'] = plot_df['Day'].astype(str).fillna('Unknown') if 'Day' in plot_df.columns else 'Unknown'
plot_df['Disease_Day'] = plot_df['Disease_Group'].astype(str) + " - " + plot_df['Day_Clean'].astype(str)

def get_day_sort_key(day_str):
    match = re.search(r'\d+', str(day_str))
    return int(match.group()) if match else 9999

plot_df['Day_Sort'] = plot_df['Day_Clean'].apply(get_day_sort_key)

# ---> STRICT PARSE BIOSCIENCES INTERSECTION LOGIC <---
print("Linking Enriched cells to surviving WT host counterparts...")
plot_df['Core_Barcode'] = adata.obs['bc_wells'].astype(str).values
plot_df['Well_Barcode_Key'] = plot_df['Day_Clean'].astype(str) + "_" + plot_df['Disease_Group'].astype(str) + "_" + plot_df['Infection_Clean'].astype(str) + "_" + plot_df['Core_Barcode'].astype(str)

wt_keys = set(plot_df[plot_df['dataset'] == 'WT']['Well_Barcode_Key'])
plot_df['Valid_Enriched'] = (plot_df['dataset'] == 'Enriched') & plot_df['Well_Barcode_Key'].isin(wt_keys)

total_enr = (plot_df['dataset'] == 'Enriched').sum()
valid_enr = plot_df['Valid_Enriched'].sum()
print(f" -> Retained {valid_enr:,} / {total_enr:,} Enriched cells for downstream Deep Dive Analytics ({(valid_enr/max(1, total_enr)*100):.1f}%)")

plot_df_umap = plot_df.copy()

if pd.api.types.is_categorical_dtype(adata.obs['majority_voting']):
    original_cell_types = list(adata.obs['majority_voting'].cat.categories)
else:
    original_cell_types = adata.obs['majority_voting'].dropna().unique().tolist()

plot_df_umap['majority_voting'] = pd.Categorical(plot_df_umap['majority_voting'], categories=original_cell_types, ordered=True)
plot_df_umap = plot_df_umap.sort_values('majority_voting')
plot_df = plot_df.sort_values(['Disease_Group', 'Day_Sort'])

# Sync cleaned columns back to adata for downstream DE analysis
adata.obs['Disease_Group'] = plot_df.loc[adata.obs_names, 'Disease_Group']
adata.obs['Infection_Clean'] = plot_df.loc[adata.obs_names, 'Infection_Clean']
adata.obs['Day_Clean'] = plot_df.loc[adata.obs_names, 'Day_Clean']
adata.obs['Day_Sort'] = plot_df.loc[adata.obs_names, 'Day_Sort']

CELLTYPE_COLOR_MAP = {
    'Tcm/Naive helper T cells': '#1f77b4', 'Tem cytotoxic T cells': '#2ca02c', 'Tem/Trm cytotoxic T cells': '#2ca02c',
    'Plasmablasts': '#d62728', 'Alveolar macrophages': '#ff7f0e',
    'NK cells': '#9467bd', 'Intestinal macrophages': '#aec7e8',
    'Regulatory T cells': '#98df8a', 'Classical monocytes': '#ff9896',
    'Memory B cells': '#fdbf6f', 'Macrophages': '#c5b0d5',
    'Tem/Effector helper T cells PD1+': '#f47d4d', 'Cycling T cells': '#75a3e1',
    'Age-associated B cells': '#299d5c', 'Naive B cells': '#9e9bc1',
    'CD16- NK cells': '#f1a243', 'Migratory DCs': '#7196aa',
    'Plasma cells': '#8bcebb', 'Tem/Temra cytotoxic T cells': '#e587bd',
    'MAIT cells': '#afd85c', 'CD16+ NK cells': '#fadd4d',
    'ILC3': '#1e9a79', 'Non-classical monocytes': '#c95b1a',
    'pDC': '#786db1', 'Unknown': '#e3e3e3'
}

b_mask = plot_df['majority_voting'].astype(str).str.contains(r'(?i)b cell|plasma|plasmablast', regex=True, na=False)
t_mask = plot_df['majority_voting'].astype(str).str.contains(r'(?i)t cell|tcm|tem|treg|regulatory|mait|cycling t', regex=True, na=False)
m_mask = plot_df['majority_voting'].astype(str).str.contains(r'(?i)monocyte|macrophage|dc\b', regex=True, na=False)

b_df, t_df, m_df = plot_df[b_mask].copy(), plot_df[t_mask].copy(), plot_df[m_mask].copy()
b_df_umap = plot_df_umap[plot_df_umap['majority_voting'].astype(str).str.contains(r'(?i)b cell|plasma|plasmablast', regex=True, na=False)].copy()
t_df_umap = plot_df_umap[plot_df_umap['majority_voting'].astype(str).str.contains(r'(?i)t cell|tcm|tem|treg|regulatory|mait|cycling t', regex=True, na=False)].copy()
m_df_umap = plot_df_umap[plot_df_umap['majority_voting'].astype(str).str.contains(r'(?i)monocyte|macrophage|dc\b', regex=True, na=False)].copy()

# Global Target Cell Map for Deep Dives & DE (Short Name -> {display_name, [categories]})
TARGET_CELLS = {
    'CD16pos_NK': {'name': 'CD16+ NK cells', 'types': ['CD16+ NK cells']},
    'ABC': {'name': 'Age-associated B cells', 'types': ['Age-associated B cells']},
    'NonClassicalMono': {'name': 'Non-classical monocytes', 'types': ['Non-classical monocytes']},
    'Naive_Helper_T': {'name': 'Tcm/Naive helper T cells', 'types': ['Tcm/Naive helper T cells']},
    'Reg_T': {'name': 'Regulatory T cells', 'types': ['Regulatory T cells']},
    'CD8_T': {'name': 'CD8+ T cells', 'types': ['Tem/Temra cytotoxic T cells', 'Tem/Trm cytotoxic T cells', 'Tem cytotoxic T cells']}
}

# ==========================================
# 3. HTML Export Helper Function
# ==========================================
def calc_box_stats(df, value_col, group_cols, compare_col='Infection_Clean'):
    results = []
    if df.empty or compare_col not in df.columns: return pd.DataFrame()
    groups = df[compare_col].unique()
    if len(groups) != 2: return pd.DataFrame()
    g1, g2 = groups[0], groups[1]
    for name, sub_df in df.groupby(group_cols):
        val1 = sub_df[sub_df[compare_col] == g1][value_col].dropna().values
        val2 = sub_df[sub_df[compare_col] == g2][value_col].dropna().values
        row = dict(zip(group_cols, name if isinstance(name, tuple) else (name,)))
        row[f'n_{g1}'] = len(val1); row[f'n_{g2}'] = len(val2)
        row[f'Mean_{g1}'] = np.mean(val1) if len(val1) > 0 else 0
        row[f'Mean_{g2}'] = np.mean(val2) if len(val2) > 0 else 0
        if len(val1) > 2 and len(val2) > 2:
            try: _, pval = stats.mannwhitneyu(val1, val2, alternative='two-sided')
            except: pval = np.nan
        else: pval = np.nan
        row['P_Value (MWU)'] = pval
        results.append(row)
    return pd.DataFrame(results).sort_values('P_Value (MWU)')

def export_plot_and_table(fig, table_df, directory, filename, title, is_umap=True, plot_height=700, bottom_margin=None):
    b_marg = bottom_margin if bottom_margin is not None else (20 if is_umap else 160)
    fig.update_layout(height=plot_height, autosize=True, plot_bgcolor='white', margin=dict(l=20, r=20, t=50, b=b_marg))
    if is_umap:
        fig.update_xaxes(showticklabels=False, title='', showgrid=False, zeroline=False)
        fig.update_yaxes(showticklabels=False, title='', showgrid=False, zeroline=False, scaleanchor="x", scaleratio=1)
    fig_html = fig.to_html(full_html=False, include_plotlyjs='cdn', default_width='100%')
    
    table_html = ""
    if table_df is not None:
        table_html = f"""
        <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 40px; border-bottom: 2px solid #dee2e6; padding-bottom: 10px;">
            <h3 style="margin: 0; border: none; padding: 0;">Summary Data</h3>
            <button class="btn-download" onclick="downloadCSV()">⬇ Download Table (CSV)</button>
        </div>
        <div class='table-container'>{table_df.to_html(classes='data-table', index=False, float_format='%.4f')}</div>
        <script>
            function downloadCSV() {{
                var table = document.querySelector(".data-table");
                if (!table) return;
                var rows = table.querySelectorAll("tr");
                var csv = [];
                for (var i = 0; i < rows.length; i++) {{
                    var row = [], cols = rows[i].querySelectorAll("td, th");
                    for (var j = 0; j < cols.length; j++) {{
                        var data = cols[j].innerText.replace(/(\\r\\n|\\n|\\r)/gm, '').replace(/"/g, '""');
                        row.push('"' + data + '"');
                    }}
                    csv.push(row.join(","));
                }}
                var csvFile = new Blob([csv.join("\\n")], {{type: "text/csv"}});
                var dl = document.createElement("a");
                dl.download = "{filename.replace('.html', '.csv')}";
                dl.href = window.URL.createObjectURL(csvFile);
                dl.style.display = "none";
                document.body.appendChild(dl);
                dl.click();
            }}
        </script>
        """
    
    html_template = f"""
    <!DOCTYPE html>
    <html><head><meta charset="utf-8"><title>{title}</title>
    <style>
        body {{ font-family: -apple-system, sans-serif; margin: 0; padding: 20px; background-color: #f8f9fa; color: #333; overflow-x: hidden; }}
        h2 {{ text-align: center; color: #2c3e50; font-weight: 600; margin-bottom: 20px; }}
        .container {{ width: 98%; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); box-sizing: border-box; }}
        .plot-container {{ display: block; width: 100%; margin-bottom: 20px; background: #fafafa; }}
        .plotly-graph-div {{ margin: 0 auto; }}
        .table-container {{ max-height: 500px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 4px; }}
        .data-table {{ width: 100%; border-collapse: collapse; font-size: 14px; }}
        .data-table th, .data-table td {{ padding: 12px 15px; border-bottom: 1px solid #dee2e6; text-align: left; }}
        .data-table th {{ background-color: #e9ecef; position: sticky; top: 0; z-index: 10; }}
        .data-table tr:hover {{ background-color: #f8f9fa; }}
        .btn-download {{ background-color: #0d6efd; color: white; border: none; padding: 8px 16px; border-radius: 4px; cursor: pointer; font-size: 14px; font-weight: 600; transition: 0.2s; }}
        .btn-download:hover {{ background-color: #0b5ed7; }}
    </style></head>
    <body><div class="container"><h2>{title}</h2><div class="plot-container">{fig_html}</div>{table_html}</div></body></html>
    """
    with open(os.path.join(directory, filename), 'w', encoding='utf-8') as f: f.write(html_template)

def sum_cat(df, col): return df.groupby(['dataset'] if col == 'dataset' else ['dataset', col]).size().reset_index(name='Total_Cells')
def sum_cont(df, col): return df.groupby(['dataset', 'majority_voting'])[col].agg(['count', 'mean', 'median', 'max']).reset_index().rename(columns={'count': 'Total_Cells'})
def sum_time(df, col): return df.groupby(['dataset', 'Disease_Group', 'Infection_Clean', 'Day_Clean', 'majority_voting'])[col].mean().reset_index().rename(columns={col: f'Mean_{col}'})

# ==========================================
# 4. Generate GLOBAL UMAPs + Tables
# ==========================================
hover_cols_global = ['majority_voting', 'Disease_Group', 'Day_Clean', 'Infection_Clean']
if 'conf_score' in plot_df.columns: hover_cols_global += ['conf_score']

if should_run(4, "global"):
    print("\n--- Generating Global Overviews ---")

    fig_ds = px.scatter(plot_df_umap, x='UMAP1', y='UMAP2', color='dataset', render_mode='webgl', title="UMAP: WT vs Enriched")
    export_plot_and_table(fig_ds, sum_cat(plot_df_umap, 'dataset'), DIR_GLOBAL, "01_Dataset.html", "Dataset Breakdown")

    # --- APPLY TARGETED RENAMING FOR GLOBAL OVERVIEW ONLY ---
    global_plot_df = plot_df_umap.copy()
    global_rename_map = {
        'Alveolar macrophages': 'Monocytes (macrophage-like)',
        'Intestinal macrophages': 'Activated / phagocytic monocyte',
        'Macrophages': 'Activated / phagocytic monocyte',
        'Non-classical monocytes': 'APC-like monocytes',
        'Tem cytotoxic T cells': 'CD8+ cytotoxic T cells',
        'CD8 cytotoxic T cells': 'CD8+ cytotoxic T cells',
        'Tem/Temra cytotoxic T cells': 'CD8+ cytotoxic T cells',
        'Tem/Trm cytotoxic T cells': 'CD8+ cytotoxic T cells'
    }
    global_plot_df['majority_voting'] = global_plot_df['majority_voting'].astype(str).replace(global_rename_map)

    # Update color map for renamed categories
    global_color_map = CELLTYPE_COLOR_MAP.copy()
    for old_name, new_name in global_rename_map.items():
        if old_name in global_color_map:
            global_color_map[new_name] = global_color_map[old_name]

    # Preserve ordering after merging categories
    global_categories = []
    for ct in original_cell_types:
        new_ct = global_rename_map.get(ct, ct)
        if new_ct not in global_categories:
            global_categories.append(new_ct)

    # 1. Standard Cell Type Distribution
    fig_ct = px.scatter(
        global_plot_df, x='UMAP1', y='UMAP2', color='majority_voting', facet_col='dataset', 
        render_mode='webgl', title="UMAP: Cell Types (WT vs Enriched)",
        color_discrete_map=global_color_map,
        category_orders={'majority_voting': global_categories},
        hover_data=hover_cols_global + (['conf_score'] if 'conf_score' in global_plot_df.columns else [])
    )
    export_plot_and_table(fig_ct, sum_cat(global_plot_df, 'majority_voting'), DIR_GLOBAL, "02_CellTypes.html", "Cell Types by Dataset")

    # ---> UPDATED: Global Stacked Bar Chart for Cell Composition (Absolute & Percentage) <---
    def generate_composition_bars(df, filename_abs, filename_pct, title_suffix=""):
        if df.empty: return
        
        wt_df = df[df['dataset'] == 'WT'].copy()
        if wt_df.empty: return
        
        def get_bar_sort(row):
            d_val = 1 if 'HC' in row['Disease_Group'] else (2 if 'Active' in row['Disease_Group'] else 3)
            t_val = int(re.search(r'\d+', row['Day_Clean']).group()) if re.search(r'\d+', row['Day_Clean']) else 99
            i_val = 1 if 'Mock' in row['Infection_Clean'] else 2
            return (d_val, t_val, i_val)
            
        wt_df['Sort_Key'] = wt_df.apply(get_bar_sort, axis=1)
        wt_df['X_Axis'] = wt_df['Disease_Group'].astype(str) + " " + wt_df['Day_Clean'].astype(str) + " " + wt_df['Infection_Clean'].astype(str)
        
        group_totals = wt_df.groupby('X_Axis').size().to_dict()
        comp_counts = wt_df.groupby(['Sort_Key', 'X_Axis', 'majority_voting']).size().reset_index(name='Count')
        comp_counts = comp_counts[comp_counts['Count'] > 0].sort_values(['Sort_Key'])
        
        # 1. Absolute Chart
        fig_abs = px.bar(
            comp_counts, x='X_Axis', y='Count', color='majority_voting', 
            title=f"Global Cell Composition (Absolute WT Counts) {title_suffix}", 
            color_discrete_map=CELLTYPE_COLOR_MAP, 
            category_orders={'majority_voting': original_cell_types, 'X_Axis': comp_counts['X_Axis'].unique()}
        )
        fig_abs.update_layout(barmode='stack')
        for x_val in comp_counts['X_Axis'].unique():
            fig_abs.add_annotation(x=x_val, y=group_totals[x_val], text=f"n={group_totals[x_val]:,}", showarrow=False, yshift=15, textangle=-45)
        fig_abs.update_xaxes(tickangle=-45, title="")
        export_plot_and_table(fig_abs, comp_counts.drop(columns=['Sort_Key']), DIR_GLOBAL, filename_abs, f"Global Composition Bar (Absolute) {title_suffix}", is_umap=False, bottom_margin=150)

        # 2. Percentage Chart
        comp_counts['Percentage'] = comp_counts.groupby('X_Axis')['Count'].transform(lambda x: x / x.sum() * 100)
        fig_pct = px.bar(
            comp_counts, x='X_Axis', y='Percentage', color='majority_voting', 
            title=f"Global Cell Composition (Relative % WT) {title_suffix}", 
            color_discrete_map=CELLTYPE_COLOR_MAP, 
            category_orders={'majority_voting': original_cell_types, 'X_Axis': comp_counts['X_Axis'].unique()}
        )
        fig_pct.update_layout(barmode='stack')
        fig_pct.update_xaxes(tickangle=-45, title="")
        export_plot_and_table(fig_pct, comp_counts.drop(columns=['Sort_Key']), DIR_GLOBAL, filename_pct, f"Global Composition Bar (Percentage) {title_suffix}", is_umap=False, bottom_margin=150)

    # Standard Composition Bars
    generate_composition_bars(plot_df, "03a_Global_Composition_Bar_Abs.html", "03b_Global_Composition_Bar_Pct.html")

    # Confidence Filtered Composition Bars
    if 'conf_score' in plot_df.columns:
        thresholds = [0.7, 0.5, 0.2]
        file_prefixes = [("03c", "03d"), ("03e", "03f"), ("03g", "03h")]
        for thresh, prefixes in zip(thresholds, file_prefixes):
            mask_conf = plot_df['conf_score'] >= thresh
            generate_composition_bars(plot_df[mask_conf], 
                                     f"{prefixes[0]}_Global_Composition_Bar_Abs_Conf_{str(thresh).replace('.', '')}.html", 
                                     f"{prefixes[1]}_Global_Composition_Bar_Pct_Conf_{str(thresh).replace('.', '')}.html", 
                                     title_suffix=f"(Conf >= {thresh})")
    else:
        print("Warning: 'conf_score' not found in data, skipping confidence-filtered composition bars.")

    fig_v_global_log = px.scatter(plot_df_umap, x='UMAP1', y='UMAP2', color='Plot_Viral_Load_Log', facet_col='dataset', hover_data=hover_cols_global, color_continuous_scale='Reds', render_mode='webgl', title="UMAP: True Viral Load (All Cells, Log1p)")

    fig_v_global_log = px.scatter(plot_df_umap, x='UMAP1', y='UMAP2', color='Plot_Viral_Load_Log', facet_col='dataset', hover_data=hover_cols_global, color_continuous_scale='Reds', render_mode='webgl', title="UMAP: True Viral Load (All Cells, Log1p)")
    export_plot_and_table(fig_v_global_log, sum_cont(plot_df_umap, 'Plot_Viral_Load_Log'), DIR_GLOBAL, "04_Global_ViralLoad_Log.html", "Global Viral Load (Log1p)")

# ==========================================
# 5. COMPARTMENT SPECIFIC UMAPS
# ==========================================
if should_run(5, "compartments"):
    def export_compartment(sub_df, sub_df_umap, prefix, name, score_col, score_name, color_scale):
        if sub_df.empty: return
        print(f"\n--- Generating {name} Dynamics ({prefix}) ---")
        fig_subtypes = px.scatter(
            sub_df_umap, x='UMAP1', y='UMAP2', color='majority_voting', facet_col='dataset', 
            hover_data=hover_cols_global, render_mode='webgl', title=f"UMAP: Subtypes ({name})",
            color_discrete_map=CELLTYPE_COLOR_MAP
        )
        export_plot_and_table(fig_subtypes, sum_cat(sub_df_umap, 'majority_voting'), DIR_COMPARTMENTS, f"{prefix}_00_Subtypes.html", f"{name} Subtypes")
        
        fig_score = px.scatter(sub_df_umap, x='UMAP1', y='UMAP2', color=score_col, facet_col='dataset', hover_data=hover_cols_global, color_continuous_scale=color_scale, render_mode='webgl', title=f"UMAP: {score_name} ({name})")
        export_plot_and_table(fig_score, sum_cont(sub_df_umap, score_col), DIR_COMPARTMENTS, f"{prefix}_01_TargetScore.html", f"{name} {score_name}")

        fig_ifn = px.scatter(sub_df_umap, x='UMAP1', y='UMAP2', color='IFN_Score', facet_col='dataset', hover_data=hover_cols_global, color_continuous_scale='Purples', render_mode='webgl', title=f"UMAP: IFN Score ({name})")
        export_plot_and_table(fig_ifn, sum_cont(sub_df_umap, 'IFN_Score'), DIR_COMPARTMENTS, f"{prefix}_02_IFN_Score.html", f"{name} IFN Response")

        dynamic_height = 300 * len(sub_df['majority_voting'].unique())
        fig_box = px.box(sub_df, x='Day_Clean', y=score_col, color='Infection_Clean', facet_col='Disease_Group', facet_row='majority_voting', points='outliers', color_discrete_map={'Mock': '#1f77b4', 'EBV': '#d62728'}, category_orders={'Disease_Group': ['HC', 'MS Active', 'MS Stable']})
        export_plot_and_table(fig_box, sum_time(sub_df, score_col), DIR_COMPARTMENTS, f"{prefix}_03_Timecourse_Score.html", f"{name} Time-Course: {score_name}", is_umap=False, plot_height=dynamic_height)

        # Targeted Contrast: HC vs Active MS Timecourse
        if prefix == "Tcell":
            print(" -> Generating Targeted T-cell Contrast (HC vs Active MS)...")
            target_mask = sub_df['Disease_Group'].isin(['HC', 'MS Active'])
            target_df = sub_df[target_mask].copy()
            
            # Combine Day and Infection for a unified X-axis
            target_df['Timecourse_Group'] = target_df['Infection_Clean'].astype(str) + " " + target_df['Day_Clean'].astype(str)
            
            # Define exact order
            tc_order = ["Mock Day 1", "EBV Day 1", "EBV Day 7", "EBV Day 15"]
            target_df = target_df[target_df['Timecourse_Group'].isin(tc_order)]
            
            # Stats calculation
            stats_list = []
            for ct in target_df['majority_voting'].unique():
                ct_df = target_df[target_df['majority_voting'] == ct]
                for dg in ['HC', 'MS Active']:
                    dg_df = ct_df[ct_df['Disease_Group'] == dg]
                    mock_vals = dg_df[dg_df['Timecourse_Group'] == "Mock Day 1"][score_col].dropna().values
                    if len(mock_vals) < 3: continue
                    
                    for ebv_tc in ["EBV Day 1", "EBV Day 7", "EBV Day 15"]:
                        ebv_vals = dg_df[dg_df['Timecourse_Group'] == ebv_tc][score_col].dropna().values
                        if len(ebv_vals) < 3: continue
                        
                        u_stat, pval = stats.mannwhitneyu(mock_vals, ebv_vals, alternative='two-sided')
                        stats_list.append({
                            'Cell_Type': ct, 'Disease_Group': dg, 'Comparison': f"{ebv_tc} vs Mock Day 1",
                            'U_Stat': u_stat, 'P_Value': pval, 'Mean_Mock': np.mean(mock_vals), 'Mean_EBV': np.mean(ebv_vals)
                        })

            stats_df = pd.DataFrame(stats_list)
            
            fig_tc_target = px.box(
                target_df, x='Timecourse_Group', y=score_col, color='Disease_Group', facet_row='majority_voting',
                category_orders={'Timecourse_Group': tc_order, 'Disease_Group': ['HC', 'MS Active']},
                color_discrete_map={'HC': '#2ca02c', 'MS Active': '#d62728'},
                title="T-cell Exhaustion/Activation: HC vs Active MS Timecourse",
                points='outliers',
                template='plotly_white'
            )
            
            # Increase height and adjust facet labels for readability
            num_subtypes = len(target_df['majority_voting'].unique())
            fig_tc_target.update_layout(
                height=350 * num_subtypes,
                margin=dict(l=50, r=20, t=80, b=100),
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
            )
            
            # Clean up facet row titles and ensure y-axis labels are clear
            fig_tc_target.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1], font=dict(size=14)))
            fig_tc_target.update_yaxes(title="Score", matches=None)
            fig_tc_target.update_xaxes(tickangle=-45, title="")
            
            export_plot_and_table(fig_tc_target, stats_df, DIR_COMPARTMENTS, "Tcell_04_Targeted_Timecourse.html", "T-cell targeted Contrast: HC vs Active MS", is_umap=False, plot_height=350 * num_subtypes)

    export_compartment(b_df, b_df_umap, "Bcell", "B-Cells", "ABC_Score", "Atypical B-Cell Signature", "RdBu_r")
    export_compartment(t_df, t_df_umap, "Tcell", "T-Cells", "T_Exhaustion_Score", "Exhaustion/Activation", "Oranges")
    export_compartment(m_df, m_df_umap, "Myeloid", "Monocytes/Macs", "Mono_Inflammatory_Score", "Inflammatory Response", "Reds")

# ==========================================
# 6. UNIVERSAL EBV TRANSCRIPT DEEP DIVE
# ==========================================
if should_run(6, "ebv"):
    print("\n--- Generating EBV Transcript Deep Dive (Valid Intersection Only) ---")
    compartments = {
        'Bcell': (b_mask, 'ABC_Score', 'ABC-High', 'ABC-Low/Zero'),
        'Tcell': (t_mask, 'T_Exhaustion_Score', 'Exhausted', 'Not Exhausted'),
        'Myeloid': (m_mask, 'Mono_Inflammatory_Score', 'Inflamed', 'Not Inflamed')
    }

    ebv_genes_in_data = [g for g in adata.var_names if g in KNOWN_EBV_GENES or g.startswith('EBV-')]

    if ebv_genes_in_data:
        for comp_prefix, (mask, target_score_col, high_label, low_label) in compartments.items():
            valid_barcodes = plot_df.index[mask & (plot_df['Plot_Viral_Load'] >= 2) & plot_df['Valid_Enriched']]
            infected_cells = adata[valid_barcodes].copy()

            if infected_cells.n_obs > 10:
                print(f" -> Processing {comp_prefix} ({infected_cells.n_obs} infected cells)...")
                
                ebv_expr = infected_cells[:, ebv_genes_in_data].X
                ebv_expr_dense = ebv_expr.toarray() if sp.issparse(ebv_expr) else np.asarray(ebv_expr)

                dd_df = pd.DataFrame(ebv_expr_dense, columns=ebv_genes_in_data, index=infected_cells.obs_names)
                dd_df['majority_voting'] = plot_df.loc[infected_cells.obs_names, 'majority_voting'].values
                dd_df['Disease_Group'] = plot_df.loc[infected_cells.obs_names, 'Disease_Group'].values
                dd_df[target_score_col] = infected_cells.obs[target_score_col].values
                dd_df['Status'] = np.where(dd_df[target_score_col] >= 0.25, high_label, low_label)

                corr_results, top_genes_list = [], []
                
                for cell_type in dd_df['majority_voting'].unique():
                    for disease in ['HC', 'MS Active', 'MS Stable']:
                        sub_df = dd_df[(dd_df['majority_voting'] == cell_type) & (dd_df['Disease_Group'] == disease)]
                        total_inf_ct = len(sub_df)
                        
                        if total_inf_ct >= 3:
                            for gene in ebv_genes_in_data:
                                gene_expr = sub_df[gene].values
                                cells_expressing = np.sum(gene_expr > 0)
                                
                                if cells_expressing > 0:
                                    facet_label = f"{cell_type} ({disease})"
                                    top_genes_list.append({
                                        'majority_voting': cell_type, 'Disease_Group': disease, 'Facet_Name': facet_label, 'EBV_Gene': gene, 
                                        'Total_Infected_Cells': total_inf_ct, '%_Expressing': (cells_expressing / total_inf_ct) * 100,
                                        'Mean_Expression_Positive_Only': np.mean(gene_expr[gene_expr > 0])
                                    })
                                    
                                    rho, pval = stats.spearmanr(gene_expr, sub_df[target_score_col].values)
                                    if not np.isnan(rho):
                                        corr_results.append({
                                            'majority_voting': cell_type, 'Disease_Group': disease, 'Facet_Name': facet_label, 'EBV_Gene': gene, 
                                            'Spearman_Rho': rho, 'P_Value': pval, 'Total_Infected_Cells': total_inf_ct, 'Cells_Expressing_Gene': cells_expressing
                                        })

                top_genes_df = pd.DataFrame(top_genes_list)
                if not top_genes_df.empty:
                    top_genes_df = top_genes_df.sort_values(['majority_voting', 'Disease_Group', 'Mean_Expression_Positive_Only'], ascending=[True, True, False])
                    top5_df = top_genes_df.groupby(['majority_voting', 'Disease_Group']).head(5)
                    
                    num_cell_types = len(top5_df['majority_voting'].unique())
                    safe_row_spacing = min(0.15, 0.5 / (num_cell_types - 1)) if num_cell_types > 1 else 0.0
                    dynamic_height = max(500, num_cell_types * 300)

                    disease_colors = {'HC': '#2ca02c', 'MS Active': '#d62728', 'MS Stable': '#1f77b4'}
                    fig_top5 = px.bar(
                        top5_df, x='Mean_Expression_Positive_Only', y='EBV_Gene', orientation='h',
                        facet_col='Disease_Group', facet_row='majority_voting', 
                        category_orders={'Disease_Group': ['HC', 'MS Active', 'MS Stable']},
                        facet_row_spacing=safe_row_spacing, facet_col_spacing=0.18, 
                        color='Disease_Group', color_discrete_map=disease_colors,
                        hover_data=['Total_Infected_Cells', '%_Expressing'], 
                        title=f"Top 5 EBV Genes in {comp_prefix} (Aligned by Disease Status)"
                    )
                    fig_top5.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
                    fig_top5.update_yaxes(matches=None, showticklabels=True, dtick=1, tickmode='linear', categoryorder='total ascending', title="")
                    fig_top5.update_xaxes(matches=None, showticklabels=True, title="Mean Expr")
                    export_plot_and_table(fig_top5, top5_df, DIR_EBV, f"{comp_prefix}_01_Top5_Genes_by_Disease.html", f"Top 5 EBV Genes per Subtype & Disease ({comp_prefix})", is_umap=False, plot_height=dynamic_height)

                corr_df = pd.DataFrame(corr_results)
                if not corr_df.empty:
                    corr_df['-log10(P-val)'] = -np.log10(corr_df['P_Value'] + 1e-300)
                    corr_df = corr_df.sort_values(['Facet_Name', 'Spearman_Rho'], ascending=[True, False])
                    
                    # Logic to avoid overlapping labels: Rank within facets and apply a minimum distance filter
                    def filter_labels(group):
                        group = group.sort_values('Spearman_Rho', ascending=False)
                        group['Display_Label'] = ""
                        last_rho, last_logp = -999, -999
                        count = 0
                        for i, row in group.iterrows():
                            if count < 8 and row['P_Value'] < 0.05:
                                # Simple overlap check: if Rho or LogP are very close to the last labeled point, skip
                                if abs(row['Spearman_Rho'] - last_rho) > 0.05 or abs(row['-log10(P-val)'] - last_logp) > 2:
                                    group.at[i, 'Display_Label'] = row['EBV_Gene']
                                    last_rho, last_logp = row['Spearman_Rho'], row['-log10(P-val)']
                                    count += 1
                        return group
                    
                    corr_df = corr_df.groupby('Facet_Name').apply(filter_labels).reset_index(drop=True)
                    
                    fig_corr = px.scatter(
                        corr_df, x='Spearman_Rho', y='-log10(P-val)', color='majority_voting', facet_col='Disease_Group', 
                        hover_data=['EBV_Gene', 'Cells_Expressing_Gene', 'majority_voting'], text='Display_Label', 
                        title=f"Viral Drivers of {target_score_col} in {comp_prefix} (Split by Disease)",
                        category_orders={'Disease_Group': ['HC', 'MS Active', 'MS Stable']}
                    )
                    fig_corr.update_traces(textposition='top center', textfont=dict(size=10))
                    fig_corr.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
                    export_plot_and_table(fig_corr, corr_df.drop(columns=['Display_Label']), DIR_EBV, f"{comp_prefix}_02_Correlations_by_Disease.html", f"Viral Drivers of {target_score_col} by Disease", is_umap=False, plot_height=600)

            # --- ADD VIRAL LOAD PLOTS FOR THIS COMPARTMENT ---
            comp_plot_df = plot_df[mask].copy()
            if not comp_plot_df.empty:
                print(f" -> Generating Viral Load box plot for {comp_prefix}...")
                fig_vl = px.box(
                    comp_plot_df, x='majority_voting', y='Plot_Viral_Load_Log', color='Infection_Clean', 
                    facet_col='Disease_Group', points='all',
                    category_orders={'Disease_Group': ['HC', 'MS Active', 'MS Stable']},
                    color_discrete_map={'Mock': '#1f77b4', 'EBV': '#d62728'},
                    title=f"Viral Load (Log1p) in {comp_prefix} Compartment"
                )
                fig_vl.update_xaxes(tickangle=-45, title="")
                fig_vl.update_layout(height=500, margin=dict(b=150))
                export_plot_and_table(fig_vl, comp_plot_df.groupby(['majority_voting', 'Disease_Group', 'Infection_Clean'])['Plot_Viral_Load'].agg(['count', 'mean', 'max']).reset_index(), 
                                      DIR_EBV, f"{comp_prefix}_03_Viral_Load_Distribution.html", f"{comp_prefix} Viral Load Distribution", is_umap=False)

# ==========================================
# 7. DAY-BY-DAY GRANULAR CONTRASTS
# ==========================================
if should_run(7, "day"):
    TARGET_DAYS = {'Day 1': 1, 'Day 7': 7, 'Day 15': 15}

    for target_day, day_int in TARGET_DAYS.items():
        print(f"\n--- Generating Day Contrasts for {target_day} ---")
        day_mask = plot_df['Day_Sort'] == day_int
        day_df = plot_df[day_mask].copy()
        if day_df.empty: continue
        day_str = target_day.replace(" ", "")

        day_df_umap = plot_df_umap[plot_df_umap['Day_Sort'] == day_int].copy()

        # Q0: Overviews
        fig_u1 = px.scatter(day_df_umap, x='UMAP1', y='UMAP2', color='majority_voting', facet_col='dataset', render_mode='webgl', title=f"{target_day}: Cell Types", color_discrete_map=CELLTYPE_COLOR_MAP, category_orders={'majority_voting': original_cell_types})
        export_plot_and_table(fig_u1, sum_cat(day_df, 'majority_voting'), DIR_DAYS, f"Q0_{day_str}_UMAP_CellTypes.html", f"{target_day} Cell Types")

        fig_u2 = px.scatter(day_df_umap, x='UMAP1', y='UMAP2', color='Disease_Group', facet_col='dataset', render_mode='webgl', title=f"{target_day}: Disease Groups", category_orders={"Disease_Group": ["HC", "MS Active", "MS Stable"]})
        export_plot_and_table(fig_u2, sum_cat(day_df, 'Disease_Group'), DIR_DAYS, f"Q0_{day_str}_UMAP_Disease.html", f"{target_day} Disease Groups")

        fig_u3 = px.scatter(day_df_umap, x='UMAP1', y='UMAP2', color='Infection_Clean', facet_col='dataset', render_mode='webgl', title=f"{target_day}: Infection Status", color_discrete_map={'Mock': '#1f77b4', 'EBV': '#d62728'})
        export_plot_and_table(fig_u3, sum_cat(day_df, 'Infection_Clean'), DIR_DAYS, f"Q0_{day_str}_UMAP_Infection.html", f"{target_day} Infection Status")

        fig_u4 = px.scatter(day_df_umap, x='UMAP1', y='UMAP2', color='Plot_Viral_Load_Log', facet_col='dataset', color_continuous_scale='Reds', render_mode='webgl', title=f"{target_day}: Viral Load (Log1p)")
        export_plot_and_table(fig_u4, sum_cont(day_df, 'Plot_Viral_Load_Log'), DIR_DAYS, f"Q0_{day_str}_UMAP_ViralLoad.html", f"{target_day} Viral Load")

        # Q1: WT Contrasts
        wt_day = day_df[day_df['dataset'] == 'WT'].copy()
        if not wt_day.empty:
            fig_q1_vl = px.box(wt_day, x='majority_voting', y='Plot_Viral_Load', color='Infection_Clean', facet_col='Disease_Group', points='outliers', category_orders={"Disease_Group": ["HC", "MS Active", "MS Stable"]}, color_discrete_map={'Mock': '#1f77b4', 'EBV': '#d62728'}, title=f"{target_day} (WT): Viral Load by Cell Type (Mock vs EBV)")
            table_vl = wt_day.groupby(['Disease_Group', 'Infection_Clean', 'majority_voting'])['Plot_Viral_Load'].agg(['count', 'mean', 'max']).reset_index()
            export_plot_and_table(fig_q1_vl, table_vl, DIR_DAYS, f"Q1_{day_str}_WT_Viral_Load.html", f"{target_day} WT Viral Load", is_umap=False, bottom_margin=180)

            comp_df = wt_day.groupby(['Disease_Group', 'Infection_Clean', 'majority_voting']).size().reset_index(name='Count')
            comp_df['Percentage'] = comp_df.groupby(['Disease_Group', 'Infection_Clean'])['Count'].transform(lambda x: x / x.sum() * 100)
            fig_q1_comp = px.bar(comp_df, x='Infection_Clean', y='Percentage', color='majority_voting', facet_col='Disease_Group', category_orders={"Disease_Group": ["HC", "MS Active", "MS Stable"]}, color_discrete_map=CELLTYPE_COLOR_MAP, title=f"{target_day} (WT): Cell Type Shifts")
            export_plot_and_table(fig_q1_comp, comp_df, DIR_DAYS, f"Q1_{day_str}_WT_Composition_Shifts.html", f"{target_day} WT Compositions", is_umap=False)

        # Q2/Q3: Enriched EBV Patterns
        enr_inf = day_df[day_df['Valid_Enriched'] & (day_df['Plot_Viral_Load'] >= 2)].copy()
        if not enr_inf.empty and ebv_genes_in_data:
            enr_adata = adata[enr_inf.index, ebv_genes_in_data].copy()
            ebv_expr_dense = enr_adata.X.toarray() if sp.issparse(enr_adata.X) else np.asarray(enr_adata.X)
            expr_df = pd.DataFrame(ebv_expr_dense, columns=ebv_genes_in_data, index=enr_adata.obs_names)
            expr_df['Disease_Group'] = enr_inf['Disease_Group'].values
            expr_df['majority_voting'] = enr_inf['majority_voting'].values
            
            grouped = expr_df.groupby(['Disease_Group', 'majority_voting'])
            mean_df = grouped.mean().reset_index().melt(id_vars=['Disease_Group', 'majority_voting'], var_name='EBV_Gene', value_name='Mean_Expr')
            frac_df = grouped.apply(lambda x: (x > 0).mean()).reset_index().melt(id_vars=['Disease_Group', 'majority_voting'], var_name='EBV_Gene', value_name='Fraction')
            dotplot_df = pd.merge(mean_df, frac_df, on=['Disease_Group', 'majority_voting', 'EBV_Gene'])
            dotplot_df = dotplot_df[dotplot_df['Fraction'] > 0.01] 
            
            if not dotplot_df.empty:
                def order_disease(group): return 1 if 'HC' in group else (2 if 'Active' in group else 3)
                dotplot_df['D_Sort'] = dotplot_df['Disease_Group'].apply(order_disease)
                dotplot_df = dotplot_df.sort_values(['D_Sort', 'majority_voting'])
                dotplot_df['Y_Axis'] = dotplot_df['majority_voting'].astype(str) + " (" + dotplot_df['Disease_Group'].astype(str) + ")"
                
                fig_q2 = px.scatter(dotplot_df, x='EBV_Gene', y='Y_Axis', size='Fraction', color='Mean_Expr', color_continuous_scale='Reds', size_max=8, title=f"{target_day} (Enriched): EBV Pattern")
                fig_q2.update_xaxes(tickangle=-90, dtick=1, tickmode='linear')
                fig_q2.update_yaxes(dtick=1, tickmode='linear')
                export_plot_and_table(fig_q2, dotplot_df.drop(columns=['D_Sort']), DIR_DAYS, f"Q2_Q3_{day_str}_Enriched_EBV_DotPlot.html", f"{target_day} Enriched EBV Patterns", is_umap=False, plot_height=900, bottom_margin=180)

        # Q4: Predictions
        q4_comp = day_df.groupby(['dataset', 'Disease_Group', 'Infection_Clean', 'majority_voting']).size().reset_index(name='Count')
        q4_comp = q4_comp[q4_comp['Count'] > 0]
        q4_comp['Percentage'] = q4_comp.groupby(['dataset', 'Disease_Group', 'Infection_Clean'])['Count'].transform(lambda x: x / x.sum() * 100)
        
        def q4_sort(row):
            d_val = 1 if 'HC' in row['Disease_Group'] else (2 if 'Active' in row['Disease_Group'] else 3)
            i_val = 1 if 'Mock' in row['Infection_Clean'] else 2
            return (d_val, i_val)
        q4_comp['Sort_Key'] = q4_comp.apply(q4_sort, axis=1)
        q4_comp = q4_comp.sort_values(['Sort_Key'])
        
        q4_comp['X_Axis'] = q4_comp['Disease_Group'].astype(str) + " " + q4_comp['Infection_Clean'].astype(str)
        fig_q4 = px.bar(q4_comp, x='X_Axis', y='Percentage', color='majority_voting', facet_col='dataset', color_discrete_map=CELLTYPE_COLOR_MAP, title=f"{target_day}: Global Cell Type Predictions", category_orders={'X_Axis': q4_comp['X_Axis'].unique()})
        fig_q4.update_xaxes(tickangle=-45)
        export_plot_and_table(fig_q4, q4_comp.drop(columns=['Sort_Key']), DIR_DAYS, f"Q4_{day_str}_Global_CellType_Predictions.html", f"{target_day} Cell Type Predictions", is_umap=False, bottom_margin=150)

# ==========================================
# 8. TARGETED DEEP DIVES (POPULATION & EBV)
# ==========================================
if should_run(8, "targeted"):
    print("\n--- Generating Targeted Deep Dives ---")

    def sort_day_disease(row):
        d_val = 1 if 'HC' in str(row['Disease_Group']) else (2 if 'Active' in str(row['Disease_Group']) else 3)
        day_str = str(row['Day_Clean'])
        match = re.search(r'\d+', day_str)
        t_val = int(match.group()) if match else 99
        return (d_val, t_val)

    for short_name, info in TARGET_CELLS.items():
        exact_cell_types = info['types']
        display_name = info['name']
        dir_name = f"Dashboards_06_{short_name}_Focus"
        os.makedirs(dir_name, exist_ok=True)
        
        # 1) Cell Population % out of Total WT Compartment
        wt_df = plot_df[plot_df['dataset'] == 'WT'].copy()
        if not wt_df.empty:
            total_wt_counts = wt_df.groupby(['Disease_Group', 'Infection_Clean', 'Day_Clean']).size().reset_index(name='Total_WT_Cells')
            ct_counts = wt_df[wt_df['majority_voting'].isin(exact_cell_types)].groupby(['Disease_Group', 'Infection_Clean', 'Day_Clean']).size().reset_index(name='CT_Cells')
            
            pop_df = pd.merge(total_wt_counts, ct_counts, on=['Disease_Group', 'Infection_Clean', 'Day_Clean'], how='left')
            pop_df['CT_Cells'] = pop_df['CT_Cells'].fillna(0)
            pop_df['Percentage_of_PBMC'] = (pop_df['CT_Cells'] / pop_df['Total_WT_Cells']) * 100
            
            # Combine for unified X-axis and filter for the requested timecourse
            pop_df['Timecourse_Group'] = pop_df['Infection_Clean'].astype(str) + " " + pop_df['Day_Clean'].astype(str)
            tc_order = ["Mock Day 1", "EBV Day 1", "EBV Day 7", "EBV Day 15"]
            pop_df = pop_df[pop_df['Timecourse_Group'].isin(tc_order)]
            
            # Prepare stats: Compare each EBV timepoint to Mock Day 1 (using Chi-square or proportional test logic)
            # For simplicity in this dashboard, we'll provide descriptive stats and can add a simple prop test
            stats_pop = []
            for dg in ['HC', 'MS Active', 'MS Stable']:
                sub_dg = pop_df[pop_df['Disease_Group'] == dg]
                mock_row = sub_dg[sub_dg['Timecourse_Group'] == "Mock Day 1"]
                if mock_row.empty: continue
                n_mock = mock_row['CT_Cells'].values[0]
                total_mock = mock_row['Total_WT_Cells'].values[0]
                
                for ebv_tc in ["EBV Day 1", "EBV Day 7", "EBV Day 15"]:
                    ebv_row = sub_dg[sub_dg['Timecourse_Group'] == ebv_tc]
                    if ebv_row.empty: continue
                    n_ebv = ebv_row['CT_Cells'].values[0]
                    total_ebv = ebv_row['Total_WT_Cells'].values[0]
                    
                    # Proportions test
                    from statsmodels.stats.proportion import proportions_ztest
                    stat, pval = proportions_ztest([n_mock, n_ebv], [total_mock, total_ebv])
                    stats_pop.append({
                        'Disease_Group': dg, 'Comparison': f"{ebv_tc} vs Mock Day 1",
                        'P_Value': pval, 'Mock_Pct': (n_mock/total_mock)*100, 'EBV_Pct': (n_ebv/total_ebv)*100
                    })
            
            stats_pop_df = pd.DataFrame(stats_pop)

            fig_pop = px.bar(
                pop_df, x='Timecourse_Group', y='Percentage_of_PBMC', color='Disease_Group', barmode='group',
                category_orders={'Timecourse_Group': tc_order, 'Disease_Group': ['HC', 'MS Active', 'MS Stable']},
                color_discrete_map={'HC': '#2ca02c', 'MS Active': '#d62728', 'MS Stable': '#1f77b4'},
                title=f"{display_name}: Population Dynamics (% of WT PBMC)",
                template='plotly_white'
            )
            fig_pop.update_xaxes(tickangle=-45, title="")
            fig_pop.update_yaxes(title="% of PBMC")
            
            # Annotate P-values on the plot (optional, might be crowded, so we'll rely on the table mostly)
            # But let's add them as text labels if requested
            export_plot_and_table(fig_pop, stats_pop_df, dir_name, f"01_{short_name}_Population_Dynamics.html", f"{display_name} Population Dynamics", is_umap=False, bottom_margin=150)

        # 2) EBV Gene Expression (Enriched dot plot + Full Raw Export)
        ct_df = plot_df[plot_df['majority_voting'].isin(exact_cell_types)].copy()
        
        ct_infected_strict = ct_df[(ct_df['Plot_Viral_Load'] >= 2) & ct_df['Valid_Enriched']]
        ct_infected_relaxed = ct_df[(ct_df['Plot_Viral_Load'] > 0) & ct_df['Valid_Enriched']]
        
        title_suffix = ""
        
        if not ct_infected_strict.empty:
            ct_infected = ct_infected_strict
            title_suffix = "(UMI ≥ 2)"
        elif not ct_infected_relaxed.empty:
            ct_infected = ct_infected_relaxed
            title_suffix = "(UMI > 0 Results)"
            print(f"   -> ⚠️ Fallback: Using UMI > 0 results for {display_name} (No UMI ≥ 2 found)")
        else:
            ct_infected = pd.DataFrame() 
        
        if not ct_infected.empty and ebv_genes_in_data:
            ct_adata = adata[ct_infected.index, ebv_genes_in_data].copy()
            ct_expr = ct_adata.X.toarray() if sp.issparse(ct_adata.X) else np.asarray(ct_adata.X)
            ct_expr_df = pd.DataFrame(ct_expr, columns=ebv_genes_in_data, index=ct_adata.obs_names)
            ct_expr_df['Disease_Group'] = ct_infected['Disease_Group'].values
            ct_expr_df['Day_Clean'] = ct_infected['Day_Clean'].values
            
            ct_grp = ct_expr_df.groupby(['Disease_Group', 'Day_Clean'])
            ct_mean = ct_grp.mean().reset_index().melt(id_vars=['Disease_Group', 'Day_Clean'], var_name='EBV_Gene', value_name='Mean_Expr')
            ct_frac = ct_grp.apply(lambda x: (x > 0).mean()).reset_index().melt(id_vars=['Disease_Group', 'Day_Clean'], var_name='EBV_Gene', value_name='Fraction')
            
            full_expr_export = pd.merge(ct_mean, ct_frac, on=['Disease_Group', 'Day_Clean', 'EBV_Gene'])
            full_expr_export = full_expr_export.rename(columns={'Mean_Expr': 'Mean_Expression', 'Fraction': 'Fraction_of_Cells_Expressing'})
            raw_csv_path = os.path.join(dir_name, f"03_{short_name}_Full_EBV_Expression_All_Genes.csv")
            full_expr_export.to_csv(raw_csv_path, index=False)

            ct_dot = full_expr_export[full_expr_export['Fraction_of_Cells_Expressing'] > 0.05].copy()
            
            if not ct_dot.empty:
                ct_dot['X_Axis'] = ct_dot['Disease_Group'].astype(str) + " (" + ct_dot['Day_Clean'].astype(str) + ")"
                ct_dot['Sort_Key'] = ct_dot.apply(sort_day_disease, axis=1)
                ct_dot = ct_dot.sort_values('Sort_Key')
                
                fig_ct_dot = px.scatter(ct_dot, x='X_Axis', y='EBV_Gene', size='Fraction_of_Cells_Expressing', color='Mean_Expression', color_continuous_scale='Reds', size_max=8, title=f"EBV Genes Expressed in {display_name} {title_suffix}")
                gene_order = ct_dot.groupby('EBV_Gene')['Mean_Expression'].sum().sort_values(ascending=True).index
                
                fig_ct_dot.update_yaxes(dtick=1, tickmode='linear', title="", categoryorder='array', categoryarray=gene_order) 
                fig_ct_dot.update_xaxes(tickangle=-45, dtick=1, tickmode='linear', categoryorder='array', categoryarray=ct_dot['X_Axis'].unique(), title="")
                
                num_genes = len(ct_dot['EBV_Gene'].unique())
                dynamic_dot_height = max(700, num_genes * 25)
                
                export_plot_and_table(fig_ct_dot, ct_dot.drop(columns=['Sort_Key']), dir_name, f"02_{short_name}_EBV_Gene_Expression.html", f"{display_name} EBV Profile", is_umap=False, plot_height=dynamic_dot_height, bottom_margin=150)

            # --- SPECIAL TARGETED ADDITIONS (Viral Load Distributions) ---
            if short_name in ["ABC", "Naive_Helper_T", "Reg_T", "CD8_T"]:
                print(f" -> Adding special {short_name} Deep Dive pages...")
                
                # Viral Load Box/Jitter Plot
                target_all_df = plot_df[plot_df['majority_voting'].isin(exact_cell_types)].copy()
                if not target_all_df.empty:
                    # Filter for the requested timecourse groups
                    tc_groups = ["Mock Day 1", "EBV Day 1", "EBV Day 7", "EBV Day 15"]
                    target_all_df['Timecourse_Group'] = target_all_df['Infection_Clean'].astype(str) + " " + target_all_df['Day_Clean'].astype(str)
                    target_all_df = target_all_df[target_all_df['Timecourse_Group'].isin(tc_groups)]
                    
                    fig_jitter = px.box(
                        target_all_df, x='Timecourse_Group', y='Plot_Viral_Load', color='Infection_Clean',
                        facet_col='Disease_Group', 
                        category_orders={'Disease_Group': ['HC', 'MS Active', 'MS Stable'], 'Timecourse_Group': tc_groups},
                        color_discrete_map={'Mock': '#1f77b4', 'EBV': '#d62728'},
                        title=f"{display_name}: Viral Load Distribution",
                        labels={'Plot_Viral_Load': 'EBV UMI Count', 'Timecourse_Group': 'Group'},
                        points='all',
                        template='plotly_white'
                    )
                    # Use matches=None to allow independent X-axes and remove empty categories per facet
                    fig_jitter.update_xaxes(matches=None, tickangle=-45)
                    export_plot_and_table(fig_jitter, target_all_df.groupby(['Disease_Group', 'Day_Clean', 'Infection_Clean'])['Plot_Viral_Load'].agg(['count', 'mean', 'max']).reset_index(), 
                                          dir_name, f"04_{short_name}_Viral_Load_Jitter.html", f"{display_name} Viral Load Distribution", is_umap=False)


                if short_name == "ABC":
                    # (ABC Only) EBV Dotplot with ALL genes
                    if ebv_genes_in_data and not ct_infected.empty:
                        ct_dot_full = full_expr_export.copy()
                        ct_dot_full['X_Axis'] = ct_dot_full['Disease_Group'].astype(str) + " (" + ct_dot_full['Day_Clean'].astype(str) + ")"
                        ct_dot_full['Sort_Key'] = ct_dot_full.apply(sort_day_disease, axis=1)
                        ct_dot_full = ct_dot_full.sort_values('Sort_Key')
                        
                        # Fill NaNs for Plotly sizing
                        ct_dot_full['Fraction_of_Cells_Expressing'] = ct_dot_full['Fraction_of_Cells_Expressing'].fillna(0.0)
                        ct_dot_full['Mean_Expression'] = ct_dot_full['Mean_Expression'].fillna(0.0)
                        
                        fig_abc_full_dot = px.scatter(
                            ct_dot_full, x='X_Axis', y='EBV_Gene', size='Fraction_of_Cells_Expressing', color='Mean_Expression', 
                            color_continuous_scale='Reds', size_max=8, title=f"ABC: Complete EBV Gene Profile (All Genes)"
                        )
                        fig_abc_full_dot.update_yaxes(dtick=1, tickmode='linear', title="")
                        fig_abc_full_dot.update_xaxes(tickangle=-45, title="")
                        export_plot_and_table(fig_abc_full_dot, ct_dot_full.drop(columns=['Sort_Key']), dir_name, "05_ABC_Full_EBV_Dotplot.html", "ABC Full EBV Profile", is_umap=False, plot_height=max(800, len(ebv_genes_in_data) * 20))

                    # (ABC Only) Raw EBV Read Counts Summary
                    if not ct_infected.empty:
                        raw_summary = ct_infected[['Disease_Group', 'Day_Clean', 'Plot_Viral_Load']].copy()
                        raw_summary = raw_summary.sort_values('Plot_Viral_Load', ascending=False)
                        fig_dummy = px.histogram(raw_summary, x='Plot_Viral_Load', nbins=50, title="ABC: Distribution of EBV UMI Counts")
                        export_plot_and_table(fig_dummy, raw_summary, dir_name, "06_ABC_Raw_EBV_Read_Counts.html", "ABC Raw EBV Read Counts", is_umap=False)


        else:
            print(f"   -> ⚠️ No valid EBV+ cells found for {display_name} (even at UMI > 0). Skipping EBV profile plots.")


# ==========================================
# 9. TARGETED DIFFERENTIAL EXPRESSION (WT ONLY)
# ==========================================
if should_run(9, "de"):
    print("\n--- Generating Targeted Differential Expression (WT base) ---")

    def run_and_plot_de(adata_obj, group_col, test_group, ref_group, filename_prefix, title):
        if sum(adata_obj.obs[group_col] == test_group) > 10 and sum(adata_obj.obs[group_col] == ref_group) > 10:
            sc.tl.rank_genes_groups(adata_obj, groupby=group_col, groups=[test_group], reference=ref_group, method='wilcoxon')
            result = adata_obj.uns['rank_genes_groups']
            
            de_df = pd.DataFrame({
                'Gene': result['names'][test_group], 'Log2FC': result['logfoldchanges'][test_group],
                'P-val': result['pvals'][test_group], 'Adj_P-val': result['pvals_adj'][test_group]
            })
            
            de_df['-log10(Adj_P-val)'] = -np.log10(de_df['Adj_P-val'] + 1e-300)
            de_df['Significance'] = 'Not Sig'
            de_df.loc[(de_df['Adj_P-val'] < 0.05) & (de_df['Log2FC'] > 1), 'Significance'] = f'Up in {test_group}'
            de_df.loc[(de_df['Adj_P-val'] < 0.05) & (de_df['Log2FC'] < -1), 'Significance'] = f'Up in {ref_group}'
            
            fig_volc = px.scatter(
                de_df, x='Log2FC', y='-log10(Adj_P-val)', color='Significance', hover_data=['Gene', 'Adj_P-val'],
                color_discrete_map={'Not Sig': '#e0e0e0', f'Up in {test_group}': '#d62728', f'Up in {ref_group}': '#1f77b4'},
                title=title
            )
            
            top_genes = de_df[de_df['Significance'] != 'Not Sig'].sort_values('Adj_P-val').head(15)
            for _, row in top_genes.iterrows():
                fig_volc.add_annotation(x=row['Log2FC'], y=row['-log10(Adj_P-val)'], text=row['Gene'], showarrow=False, yshift=10)
                
            export_plot_and_table(fig_volc, de_df.head(300), DIR_DE, f"{filename_prefix}_Volcano.html", title, is_umap=False)
            
            up_genes = de_df[de_df['Significance'] == f'Up in {test_group}'].sort_values('Adj_P-val')['Gene'].tolist()
            down_genes = de_df[de_df['Significance'] == f'Up in {ref_group}'].sort_values('Adj_P-val')['Gene'].tolist()
            
            enrich_df = pd.DataFrame({f'Up_in_{test_group}': pd.Series(up_genes), f'Up_in_{ref_group}': pd.Series(down_genes)})
            enrich_df.to_csv(os.path.join(DIR_DE, f"{filename_prefix}_Functional_Enrichment_Genes.csv"), index=False)

    for short_name, info in TARGET_CELLS.items():
        exact_cell_types = info['types']
        display_name = info['name']
        print(f" -> Running strict WT contrasts for {display_name}...")
        mask = (adata.obs['majority_voting'].isin(exact_cell_types)) & (adata.obs['dataset'] == 'WT')
        adata_sub = adata[mask].copy()
        
        if adata_sub.n_obs < 20: continue

        adata_sub.obs['DE_Combo'] = adata_sub.obs['Disease_Group'].astype(str) + "_" + adata_sub.obs['Infection_Clean'].astype(str) + "_" + adata_sub.obs['Day_Clean'].astype(str)

        run_and_plot_de(adata_sub, 'DE_Combo', 'MS Active_Mock_Day 1', 'HC_Mock_Day 1', f"01_{short_name}_MS_vs_HC_Mock_D1", f"{display_name}: Active MS vs HC (Day 1 Mock)")
        run_and_plot_de(adata_sub, 'DE_Combo', 'MS Active_EBV_Day 7', 'MS Active_EBV_Day 1', f"02_{short_name}_ActiveMS_EBV_D7_vs_D1", f"{display_name}: Active MS EBV (Day 7 vs Day 1)")

# ==========================================
# 10. PATHWAY DYNAMICS (FUNCTIONAL ENRICHMENT)
# ==========================================
if should_run(10, "pathway"):
    print("\n--- Generating Pathway Dynamics (Memory-Safe) ---")

    COMPARTMENTS_PATHWAY = {
        'B-Cells': 'b cell|b|plasma|plasmablast',
        'T-Cells': 't cell|tcm|tem|treg|regulatory|cycling t',
        'Myeloid': 'monocyte|macrophage'
    }

    TARGET_DAYS_MAP = {1: 'Day 1', 7: 'Day 7', 15: 'Day 15'}
    DISEASES = ['HC', 'MS Active', 'MS Stable']
    INFECTIONS = ['Mock', 'EBV', 'Combined']
    GENE_SETS = ['GO_Biological_Process_2023', 'KEGG_2021_Human']

    wt_mask = adata.obs['dataset'] == 'WT'
    print(f"Identified {wt_mask.sum():,} WT cells for Pathway Analysis.")

    master_results = []
    html_summary = []

    for comp_name, comp_regex in COMPARTMENTS_PATHWAY.items():
        print(f"================ Processing {comp_name} ================")
        
        comp_mask = wt_mask & adata.obs['majority_voting'].astype(str).str.contains(comp_regex, case=False, na=False)
        if comp_mask.sum() == 0: continue
            
        for day_int, day_str in TARGET_DAYS_MAP.items():
            for disease in DISEASES:
                for inf in INFECTIONS:
                    
                    strat_name = f"{day_str}-{disease}-{inf}"
                    
                    m_day = adata.obs['Day_Sort'] == day_int
                    m_dis = adata.obs['Disease_Group'] == disease
                    
                    if inf == 'Combined':
                        m_inf = adata.obs['Infection_Clean'].isin(['Mock', 'EBV'])
                    else:
                        m_inf = adata.obs['Infection_Clean'] == inf
                        
                    m_target_full = comp_mask & m_day & m_dis & m_inf
                    m_bg_full = comp_mask & ~m_target_full
                    
                    n_cells = m_target_full.sum()
                    n_bg = m_bg_full.sum()
                    
                    if n_cells < 20 or n_bg < 20:
                        continue
                        
                    print(f" -> Analyzing: {strat_name} (Target: {n_cells:,} cells | Background: {n_bg:,} cells)")
                    
                    target_indices = np.where(m_target_full)[0]
                    bg_indices = np.where(m_bg_full)[0]
                    
                    MAX_CELLS = 1000
                    if len(target_indices) > MAX_CELLS:
                        np.random.seed(42)
                        target_indices = np.random.choice(target_indices, MAX_CELLS, replace=False)
                    if len(bg_indices) > MAX_CELLS:
                        np.random.seed(42)
                        bg_indices = np.random.choice(bg_indices, MAX_CELLS, replace=False)
                        
                    subset_idx = np.concatenate([target_indices, bg_indices])
                    adata_tiny = adata[subset_idx].copy()
                    adata_tiny.obs['Comparison'] = ['Target']*len(target_indices) + ['Background']*len(bg_indices)
                    
                    try:
                        sc.tl.rank_genes_groups(adata_tiny, groupby='Comparison', groups=['Target'], reference='Background', method='t-test', use_raw=True)
                        
                        res = adata_tiny.uns['rank_genes_groups']
                        genes = res['names']['Target']
                        pvals = res['pvals_adj']['Target']
                        logfcs = res['logfoldchanges']['Target']
                        
                        sig_mask = (pvals < 0.05) & (logfcs > 0.5)
                        top_genes = list(genes[sig_mask][:150])
                        
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
                                    
                                    if idx < 3: 
                                        pathway_html_list.append(f"• <b>{row['Term_Clean']}</b><br><i>&nbsp;&nbsp;p={row['Adjusted P-value']:.2e} | Overlap: {over_genes}/{total_path_genes}</i>")
                                        matched_genes_list = str(row['Genes']).split(';')
                                        top_20 = ", ".join(matched_genes_list[:20])
                                        if len(matched_genes_list) > 20: top_20 += "..."
                                        genes_html_list.append(f"<div style='margin-bottom: 12px;'><b>{row['Term_Clean']}</b>:<br><span style='font-size: 0.85em; color: #555; line-height: 1.2;'>{top_20}</span></div>")
                                        
                                html_summary.append({
                                    'Compartment': comp_name, 'Stratification': f"<b>{strat_name}</b>",
                                    'Group Cells': f"{n_cells:,}", 'Input Marker Genes': num_input_genes,
                                    'Top Active Pathways': "<br><br>".join(pathway_html_list),
                                    'Top 20 Overlapping Genes': "".join(genes_html_list)
                                })
                            else:
                                html_summary.append({'Compartment': comp_name, 'Stratification': f"<b>{strat_name}</b>", 'Group Cells': f"{n_cells:,}", 'Input Marker Genes': num_input_genes, 'Top Active Pathways': "<i>No significant pathways identified.</i>", 'Top 20 Overlapping Genes': "<i>N/A</i>"})
                            time.sleep(1.0)
                        else:
                            html_summary.append({'Compartment': comp_name, 'Stratification': f"<b>{strat_name}</b>", 'Group Cells': f"{n_cells:,}", 'Input Marker Genes': num_input_genes, 'Top Active Pathways': "<i>Not enough marker genes (<10) for pathway analysis.</i>", 'Top 20 Overlapping Genes': "<i>N/A</i>"})
                            
                    except Exception as e:
                        print(f"      Error testing {strat_name}: {e}")

                    del adata_tiny

    df_master = pd.DataFrame(master_results)
    if not df_master.empty:
        csv_path = os.path.join(DIR_PATHWAYS, "Master_Pathway_Stratifications_WT.csv")
        df_master.to_csv(csv_path, index=False)
        print(f"\nSaved complete deep-stat data to {csv_path}")

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
        
        html_path = os.path.join(DIR_PATHWAYS, "Interactive_Pathway_Dynamics_WT.html")
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_template)
        print(f"Saved interactive dashboard to {html_path}")

print("\nMaster Dashboard Generation Complete! Check all output directories.")