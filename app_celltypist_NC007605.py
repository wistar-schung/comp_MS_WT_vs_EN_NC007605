import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import os
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import traceback
import re
import matplotlib.pyplot as plt
import scipy.sparse as sp
import mygene

st.set_page_config(page_title="EBV WT Analysis", layout="wide")
st.title("🦠 EBV WT Viral Load Dashboard")

# ==============================================================================
# 1. CONSTANTS & MARKERS
# ==============================================================================
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

MS_MARKER_LIST = ['HLA-DRB1', 'IL7R', 'C1QA', 'CR1', 'DYSF', 'ZNF638', 'MIF', 'DDT', 'EOMES', 'TBX21', 'CR2', 'CXCR3', 'CXCR4', 'FOXP3']
LATENT_GENES = ['LMP-1', 'LMP-2A', 'LMP-2B', 'EBNA-1', 'EBNA-2', 'EBNA-3A', 'EBNA-3B', 'EBNA-3C', 'EBNA-LP']
LYTIC_GENES  = ['BZLF1', 'BZLF2', 'LF1', 'LF3']

MARKER_DICT = {
    'Naive CD4+ T': ['IL7R', 'CCR7', 'TCF7'],  
    'Memory CD4+ T': ['IL7R', 'S100A4', 'CD40LG'], 
    'Tregs': ['FOXP3', 'IL2RA', 'CTLA4', 'CXCR3', 'CXCR4'], 
    'CD8+ T': ['CD8A', 'CD8B', 'GZMK'], 
    'Naive CD8+ T': ['CD8A', 'CCR7', 'LEF1'], 
    'NK Cells': ['GNLY', 'NKG7', 'KLRB1', 'FCGR3A', 'TYROBP'], 
    'Naive B': ['MS4A1', 'CD79A', 'IGHD', 'TCL1A'], 
    'Memory B': ['MS4A1', 'CD79A', 'CD27', 'AIM2'], 
    'Plasma Cells': ['MZB1', 'JCHAIN', 'IGHG1', 'XBP1'], 
    'CD14+ Monocytes': ['CD14', 'LYZ', 'VCAN'], 
    'FCGR3A+ Monocytes': ['FCGR3A', 'MS4A7', 'CD14'], 
    'mDC': ['FCER1A', 'CST3', 'CLEC10A', 'HLA-DQA1'], 
    'pDC': ['IL3RA', 'LILRA4', 'CLEC4C'], 
    'Progenitors': ['CD34', 'KIT'],
    'Platelets': ['PPBP', 'PF4']
}

@st.cache_data
def load_bcell_markers(csv_path="B-Cell_Gene_Population.csv"):
    try:
        if not os.path.exists(csv_path): return {}
        df_markers = pd.read_csv(csv_path, skiprows=2)
        expected_cols = ['cluster', 'p_val_adj', 'avg_log2FC', 'gene']
        if any(col not in df_markers.columns for col in expected_cols): return {}

        cleaned_genes = []
        seen_genes = set()
        for g in df_markers['gene'].astype(str):
            g_clean = g
            match = re.search(r'(\d+)$', g)
            if match:
                digits = match.group(1)
                for i in range(1, len(digits) + 1):
                    possible_base = g[:-i]
                    if possible_base in seen_genes:
                        g_clean = possible_base
                        break
            seen_genes.add(g_clean)
            seen_genes.add(g) 
            cleaned_genes.append(g_clean)
            
        df_markers['gene_clean'] = cleaned_genes
        df_markers = df_markers.sort_values(by=['cluster', 'p_val_adj', 'avg_log2FC'], ascending=[True, True, False])
        
        bcell_dict = {}
        for cluster in df_markers['cluster'].dropna().unique():
            top_genes = df_markers[df_markers['cluster'] == cluster]['gene_clean'].drop_duplicates().head(10).tolist()
            bcell_dict[cluster] = top_genes
        return bcell_dict
    except Exception as e:
        return {}

BCELL_MARKERS = load_bcell_markers()

def get_extended_palette():
    # --- True Publication-Ready Palette (ColorBrewer Mix) ---
    return [
        "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", 
        "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", 
        "#FD7446", "#709AE1", "#31A354", "#9E9AC8", "#F59D3D",
        "#729EAC", "#8DD3C7", "#E78AC3", "#A6D854", "#FFD92F",
        "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
        "#E5C494", "#B3B3B3", "#8DA0CB", "#66C2A5", "#FC8D62"
    ]

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

def classify_viral_state(adata):
    valid_latent = [g for g in LATENT_GENES if g in adata.var_names]
    valid_lytic  = [g for g in LYTIC_GENES if g in adata.var_names]
    
    if valid_latent: sc.tl.score_genes(adata, valid_latent, score_name='Latent_Score', use_raw=False)
    else: adata.obs['Latent_Score'] = 0.0
    if valid_lytic: sc.tl.score_genes(adata, valid_lytic, score_name='Lytic_Score', use_raw=False)
    else: adata.obs['Lytic_Score'] = 0.0

    def assign_state(row):
        if row['Viral_Counts'] < 2: return "Negative"
        elif row['Lytic_Score'] > row['Latent_Score']: return "Lytic-Like"
        else: return "Latent-Like"
    adata.obs['Viral_State'] = adata.obs.apply(assign_state, axis=1)
    return adata

def get_day_sort_key(day_str):
    try:
        match = re.search(r'\d+', str(day_str))
        if match: return int(match.group())
        return 9999
    except: return 9999

# ==============================================================================
# 2. DATA LOADING & DYNAMIC RECALCULATION
# ==============================================================================
@st.cache_resource
def load_preprocessed_data(min_genes_cutoff):
    file_path = os.path.join("optimized_data", f"data_{min_genes_cutoff}_annotated_post_qc.h5ad")
    if not os.path.exists(file_path):
        return None, f"File not found: {file_path}. Run preprocess.py first.", None

    try:
        adata = sc.read_h5ad(file_path)
        if not sp.issparse(adata.X): adata.X = sp.csr_matrix(adata.X)

        # ---> NEW: Intercept and apply PBMC Relabeling on the fly <---
        relabel_map = {
            'Alveolar macrophages': 'Alveolar macrophages',
            'Intestinal macrophages': 'Intestinal macrophages',
            'ILC3': 'Innate Lymphoid Cells',
            'Tem/Trm cytotoxic T cells': 'Tem cytotoxic T cells',
            'Tem/Temra cytotoxic T cells': 'Tem cytotoxic T cells' # Also catching Temra based on your screenshot!
        }
        
        if 'majority_voting' in adata.obs.columns:
            adata.obs['majority_voting'] = adata.obs['majority_voting'].replace(relabel_map)
        if 'Predicted_Cell_Type' in adata.obs.columns:
            adata.obs['Predicted_Cell_Type'] = adata.obs['Predicted_Cell_Type'].replace(relabel_map)
        # ------------------------------------------------------------
        
        # --- FIX: Dynamically Recalculate Viral Counts from Raw ---
        adata.var_names = adata.var_names.str.replace(r'(?i)_type1$', '', regex=True)
        adata.var_names_make_unique()
        
        ebv_genes = [g for g in adata.var_names if g in EXACT_EBV_GENES or g.startswith('EBV-')]
        
        if adata.raw is not None:
            raw_var_names = adata.raw.var_names.str.replace(r'(?i)_type1$', '', regex=True)
            ebv_raw_indices = [i for i, g in enumerate(raw_var_names) if g in ebv_genes]
            ebv_matrix = adata.raw.X[:, ebv_raw_indices]
        else:
            ebv_matrix = adata[:, ebv_genes].X
            
        # Overwrite the old, inaccurate counts with the new strict GTF counts
        adata.obs['Viral_Counts'] = np.ravel(ebv_matrix.sum(axis=1))
        adata.obs['Viral_Counts_Log'] = np.log1p(adata.obs['Viral_Counts'])
        # ------------------------------------------------------------
        
        sens_data = adata.uns.get('sensitivity_df', None)
        sensitivity_df = pd.DataFrame(sens_data) if sens_data is not None else None
            
        adata = classify_viral_state(adata)
            
        if 'Day' in adata.obs.columns:
            adata.obs['Day_Sorted'] = adata.obs['Day'].astype(str).apply(get_day_sort_key).astype(int)
            min_day = adata.obs['Day_Sorted'].min()
            max_day = adata.obs['Day_Sorted'].max()
            
            if 'dpt_pseudotime' not in adata.obs.columns:
                if 'neighbors' not in adata.uns: sc.pp.neighbors(adata)
                sc.tl.diffmap(adata)
                
                day1_mask = adata.obs['Day_Sorted'] == min_day
                day1_idx = np.flatnonzero(day1_mask)
                
                if len(day1_idx) > 0:
                    dc_centroid = adata.obsm['X_diffmap'][day1_idx, 1:4].mean(axis=0)
                    distances = np.linalg.norm(adata.obsm['X_diffmap'][day1_idx, 1:4] - dc_centroid, axis=1)
                    adata.uns['iroot'] = day1_idx[np.argmin(distances)]
                    sc.tl.dpt(adata)
            
            if 'dpt_pseudotime' in adata.obs.columns:
                median_dpt_min_day = adata.obs.loc[adata.obs['Day_Sorted'] == min_day, 'dpt_pseudotime'].median()
                median_dpt_max_day = adata.obs.loc[adata.obs['Day_Sorted'] == max_day, 'dpt_pseudotime'].median()
                
                if median_dpt_max_day < median_dpt_min_day:
                    dpt_max_val = adata.obs['dpt_pseudotime'].max()
                    adata.obs['dpt_pseudotime'] = dpt_max_val - adata.obs['dpt_pseudotime']
            
        return adata, ebv_genes, sensitivity_df

    except Exception as e:
        return None, f"{str(e)}\n\n{traceback.format_exc()}", None

# ==============================================================================
# 3. APP SETUP & SIDEBAR
# ==============================================================================
with st.sidebar:
    st.header("⚙️ Settings")
    min_genes_choice = st.radio("Min Genes Filter:", [200, 50], index=0)
    
    st.markdown("---")
    st.header("🏷️ Annotation Method")
    annotation_source = st.radio("Cell Annotation Source:", ["CellTypist (majority_voting)", "Manual Markers (Predicted_Cell_Type)"])
    
    st.markdown("---")
    res_choice = st.radio("Clustering Res:", [0.2, 0.5, 1.0], index=2, horizontal=True)
    min_cluster_size = st.number_input("Min Cluster Size:", 0, value=0, step=10)
    st.markdown("---")
    st.header("🔍 Gene Search")
    search_term = st.text_input("Gene name:", "EBNA")
    
with st.spinner(f"Loading Data (> {min_genes_choice} genes)..."):
    adata_obj, result_second, sensitivity_df = load_preprocessed_data(min_genes_choice)

if adata_obj is None:
    st.error("❌ Analysis Failed"); st.error(result_second); st.stop()

adata = adata_obj
found_viral_genes = result_second

CT_COL = 'majority_voting' if 'CellTypist' in annotation_source else 'Predicted_Cell_Type'

if CT_COL not in adata.obs.columns:
    st.sidebar.error(f"⚠️ Column '{CT_COL}' not found. Falling back to existing labels.")
    CT_COL = 'Predicted_Cell_Type' if 'Predicted_Cell_Type' in adata.obs.columns else 'leiden'

keep_cols = [
    CT_COL, 'Predicted_Cell_Type', 'majority_voting', 'predicted_labels', 'Infection', 'Viral_State', 
    'Latent_Score', 'Lytic_Score', 'leiden', 'Disease_Condition (Detail)', 'Day', 'Day_Cat', 
    'Viral_Counts', 'Viral_Counts_Log', 'Pseudotime', 'dpt_pseudotime'
]

keep_cols = list(dict.fromkeys(keep_cols))
adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]

existing_cols = [c for c in keep_cols if c in adata.obs.columns]
adata.obs = adata.obs[existing_cols].copy()

if f'leiden_res_{res_choice}' in adata.obs.columns: 
    adata.obs['leiden'] = adata.obs[f'leiden_res_{res_choice}']
    with st.spinner("Mapping PAGA connections for current clusters..."):
        sc.tl.paga(adata, groups='leiden')

st.sidebar.markdown("---")
st.sidebar.write(f"**Viral Genes Found:** {len(found_viral_genes)}")
if len(found_viral_genes) > 0:
    with st.sidebar.expander("See full list"): st.write(found_viral_genes)

st.subheader("📊 Dashboard")
n_cells = adata.shape[0]
n_infected = sum(adata.obs['Viral_Counts'] >= 2)
m1, m2, m3 = st.columns(3)
with m1: st.metric("Total Cells", f"{n_cells:,}")
with m2: st.metric("Infected Cells", f"{n_infected:,}")
with m3: st.metric("Viral Reads", f"{int(adata.obs['Viral_Counts'].sum()):,}")

# ==============================================================================
# 4. UMAP VISUALIZATION
# ==============================================================================
plot_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
safe_plot_cols = []

priority = [CT_COL, 'Infection', 'Viral_State', 'leiden', 'Disease_Condition (Detail)', 'Day', 'MS Markers', 'QC Metrics', 'Viral_Counts', 'B-Cell Analysis', 'Disease Comparison', 'Mock Breakdown', 'Trajectory Analysis']

for p in priority:
    if p in adata.obs.columns or p in ['MS Markers', 'QC Metrics', 'B-Cell Analysis', 'Disease Comparison', 'Mock Breakdown', 'Trajectory Analysis']: 
        safe_plot_cols.append(p)

for c in adata.obs.columns:
    if c not in safe_plot_cols and c not in ['n_genes', 'total_counts', 'n_genes_by_counts', 'total_counts_mt', 'pct_counts_mt'] and not c.startswith('leiden_res'):
        safe_plot_cols.append(c)

for c in safe_plot_cols:
    if c in adata.obs.columns: plot_df[c] = adata.obs[c].values

if 'Viral_Counts' in adata.obs: plot_df['Viral_Counts'] = adata.obs['Viral_Counts'].values
if 'Viral_Counts_Log' in adata.obs: plot_df['Viral_Counts_Log'] = adata.obs['Viral_Counts_Log'].values

if len(plot_df) > 50000:
    st.sidebar.warning(f"⚠️ Displaying 50,000 / {len(plot_df):,} cells for performance.")
    plot_df = plot_df.sample(50000, random_state=42)

c1, c2 = st.columns([2, 1])
with c1: color_by = st.selectbox("Color UMAP By:", safe_plot_cols, index=0)

PLOT_WIDTH = 800  
PLOT_HEIGHT = 700 

if color_by == 'Viral_Counts':
    show_log_side = st.checkbox("Show Log-Scale", value=False)
    if show_log_side:
        u1, u2 = st.columns(2)
        u1.plotly_chart(px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Viral_Counts', title="Linear", render_mode='webgl', color_continuous_scale='Reds', height=PLOT_HEIGHT), use_container_width=True)
        u2.plotly_chart(px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Viral_Counts_Log', title="Log1p", render_mode='webgl', color_continuous_scale='Reds', height=PLOT_HEIGHT), use_container_width=True)
    else:
        st.plotly_chart(px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Viral_Counts', title="Viral Counts (Linear)", render_mode='webgl', color_continuous_scale='Reds', width=PLOT_WIDTH, height=PLOT_HEIGHT), use_container_width=False)

elif color_by == 'B-Cell Analysis':
    b_mask = plot_df[CT_COL].astype(str).str.contains('b cell|b|plasma', case=False, na=False)
    plot_df['Is_B_Cell'] = np.where(b_mask, plot_df[CT_COL], 'Other')
    fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Is_B_Cell', render_mode='webgl', title="UMAP: B-Cell Populations", color_discrete_map={'Other': '#eeeeee'}, width=PLOT_WIDTH, height=PLOT_HEIGHT)
    st.plotly_chart(fig, use_container_width=False)

elif color_by == 'Disease Comparison':
    fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Disease_Condition (Detail)', render_mode='webgl', title="UMAP: Disease Conditions", width=PLOT_WIDTH, height=PLOT_HEIGHT)
    st.plotly_chart(fig, use_container_width=False)

elif color_by == 'Mock Breakdown':
    is_mock = plot_df['Infection'].astype(str).str.contains('Mock', case=False, na=False)
    plot_df['Is_Mock'] = np.where(is_mock, 'Mock', 'Other')
    fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Is_Mock', render_mode='webgl', color_discrete_map={'Other': '#eeeeee', 'Mock': '#1f77b4'}, title="UMAP: Mock Highlighted", width=PLOT_WIDTH, height=PLOT_HEIGHT)
    st.plotly_chart(fig, use_container_width=False)

elif color_by == 'Trajectory Analysis':
    pseudo_col = 'dpt_pseudotime' if 'dpt_pseudotime' in plot_df.columns else 'Pseudotime' if 'Pseudotime' in plot_df.columns else None
    
    u1, u2 = st.columns(2)
    if 'Day' in plot_df.columns:
        fig1 = px.scatter(plot_df, x='UMAP1', y='UMAP2', color='Day', render_mode='webgl', title="UMAP: Actual Timepoints", height=500)
        u1.plotly_chart(fig1, use_container_width=True)
    if pseudo_col:
        fig2 = px.scatter(plot_df, x='UMAP1', y='UMAP2', color=pseudo_col, render_mode='webgl', color_continuous_scale='Viridis', title="UMAP: Pseudotime Progression", height=500)
        u2.plotly_chart(fig2, use_container_width=True)
    elif not pseudo_col and 'Day' not in plot_df.columns:
        st.warning("Missing required Temporal/Pseudotime data.")

elif color_by not in ['MS Markers', 'QC Metrics', 'B-Cell Analysis', 'Disease Comparison', 'Mock Breakdown', 'Trajectory Analysis']:
    if pd.api.types.is_numeric_dtype(plot_df[color_by]):
        fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color=color_by, render_mode='webgl', color_continuous_scale='Reds', width=PLOT_WIDTH, height=PLOT_HEIGHT)
    else:
        # ---> APPLY THE STRICT DICTIONARY IF COLORING BY CELL TYPE <---
        if color_by == CT_COL or color_by == 'majority_voting':
            fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color=color_by, render_mode='webgl', color_discrete_map=CELLTYPE_COLOR_MAP, width=PLOT_WIDTH, height=PLOT_HEIGHT)
        else:
            # Fall back to the automatic extended palette for things like 'leiden' or 'Day'
            fig = px.scatter(plot_df, x='UMAP1', y='UMAP2', color=color_by, render_mode='webgl', color_discrete_sequence=get_extended_palette(), width=PLOT_WIDTH, height=PLOT_HEIGHT)
    st.plotly_chart(fig, use_container_width=False)

# ==============================================================================
# 5. DETAILED SECTIONS
# ==============================================================================

if color_by == 'Viral_Counts':
    st.markdown("---")
    st.subheader("🔬 Viral Gene Breakdown")
    if len(found_viral_genes) > 0:
        sel = st.multiselect("Select genes:", found_viral_genes, default=found_viral_genes[:3])
        if sel:
            sub_adata = adata[plot_df.index] 
            n_rows = (len(sel)+2)//3
            fig_sub = make_subplots(rows=n_rows, cols=3, subplot_titles=sel)
            bx = sub_adata.obsm['X_umap'][:, 0]; by = sub_adata.obsm['X_umap'][:, 1]
            for i, gene in enumerate(sel):
                expr = sub_adata[:, gene].X.toarray().flatten() if sp.issparse(sub_adata.X) else np.asarray(sub_adata[:, gene].X).flatten()
                fig_sub.add_trace(go.Scattergl(x=bx, y=by, mode='markers', marker=dict(color='#e0e0e0', size=2), showlegend=False), row=(i//3)+1, col=(i%3)+1)
                idx = expr > 0
                if sum(idx)>0: fig_sub.add_trace(go.Scattergl(x=bx[idx], y=by[idx], mode='markers', marker=dict(color=expr[idx], colorscale='Reds', size=3), showlegend=False), row=(i//3)+1, col=(i%3)+1)
            fig_sub.update_layout(height=300*n_rows)
            st.plotly_chart(fig_sub, use_container_width=True)

if color_by == 'Viral_State':
    st.markdown("---")
    st.subheader("🧬 Viral Status Analysis (Latent vs Lytic)")
    
    infected_df = plot_df[plot_df['Viral_State'] != 'Negative'].copy()
    
    if infected_df.empty:
        st.info("No infected cells found to classify.")
    else:
        t1, t2, t3 = st.tabs(["Composition by Cell Type", "Composition by Condition", "Latent/Lytic Scores"])
        
        with t1:
            fig1 = px.histogram(infected_df, x=CT_COL, color='Viral_State', barnorm='percent', 
                                title="Latent vs Lytic Proportions by Cell Type", barmode='stack',
                                color_discrete_map={'Latent-Like': '#1f77b4', 'Lytic-Like': '#d62728'})
            st.plotly_chart(fig1, use_container_width=True)
            
        with t2:
            if 'Disease_Condition (Detail)' in infected_df.columns:
                fig2 = px.histogram(infected_df, x='Disease_Condition (Detail)', color='Viral_State', barnorm='percent', 
                                    title="Latent vs Lytic Proportions by Condition", barmode='stack',
                                    color_discrete_map={'Latent-Like': '#1f77b4', 'Lytic-Like': '#d62728'})
                st.plotly_chart(fig2, use_container_width=True)
        
        with t3:
            if 'Latent_Score' in adata.obs.columns and 'Lytic_Score' in adata.obs.columns:
                score_df = adata.obs.loc[infected_df.index, [CT_COL, 'Viral_State', 'Latent_Score', 'Lytic_Score']].copy()
                
                c1, c2 = st.columns(2)
                c1.plotly_chart(px.box(score_df, x=CT_COL, y='Latent_Score', color='Viral_State', title="Latent Score Distribution", points='all'), use_container_width=True)
                c2.plotly_chart(px.box(score_df, x=CT_COL, y='Lytic_Score', color='Viral_State', title="Lytic Score Distribution", points='all'), use_container_width=True)
            else:
                st.info("Latent/Lytic scores not found.")

if color_by == 'MS Markers':
    st.markdown("---")
    valid_ms = [g for g in MS_MARKER_LIST if g in adata.var_names]
    if valid_ms:
        t1, t2 = st.tabs(["Grid View", "Dot Plot"])
        with t1:
            sub_adata = adata[plot_df.index] 
            n = len(valid_ms); cols=3; rows=(n+cols-1)//cols
            fig = make_subplots(rows=rows, cols=cols, subplot_titles=valid_ms)
            bx = sub_adata.obsm['X_umap'][:, 0]; by = sub_adata.obsm['X_umap'][:, 1]
            for i, g in enumerate(valid_ms):
                expr = sub_adata[:, g].X.toarray().flatten() if sp.issparse(sub_adata.X) else np.asarray(sub_adata[:, g].X).flatten()
                fig.add_trace(go.Scattergl(x=bx, y=by, mode='markers', marker=dict(color='#e0e0e0', size=2), showlegend=False), row=(i//cols)+1, col=(i%cols)+1)
                idx = expr > 0
                if sum(idx)>0: fig.add_trace(go.Scattergl(x=bx[idx], y=by[idx], mode='markers', marker=dict(color=expr[idx], colorscale='Viridis', size=3), showlegend=False), row=(i//cols)+1, col=(i%cols)+1)
            fig.update_layout(height=300*rows)
            st.plotly_chart(fig, use_container_width=True)
        with t2:
            dot_res = []
            for ct in sorted(plot_df[CT_COL].dropna().unique()):
                sub = adata[adata.obs[CT_COL]==ct]
                if sub.n_obs>0:
                    for g in valid_ms:
                        e = sub[:, g].X.toarray().flatten() if sp.issparse(sub.X) else np.asarray(sub[:, g].X).flatten()
                        dot_res.append({'Cell Type': ct, 'Gene': g, 'Mean': np.mean(e), 'Pct': np.mean(e>0)*100})
            st.plotly_chart(px.scatter(pd.DataFrame(dot_res), x='Gene', y='Cell Type', size='Pct', color='Mean', color_continuous_scale='Reds', render_mode='webgl'), use_container_width=True)

if color_by == 'Infection':
    st.markdown("---")
    tab1, tab2 = st.tabs(["🧬 Cell Type Composition", "🦠 Viral Load by Condition"])
    with tab1:
        st.plotly_chart(px.histogram(plot_df, x='Infection', color=CT_COL, barnorm='percent', title="Cell Type Proportions"), use_container_width=True)
    with tab2:
        df_pos = plot_df[plot_df['Viral_Counts'] >= 2]
        if not df_pos.empty:
            st.plotly_chart(px.box(df_pos, x='Disease_Condition (Detail)', y='Viral_Counts', color='Disease_Condition (Detail)', points='all', title="Viral Load (Infected Only)"), use_container_width=True)

if color_by == 'B-Cell Analysis':
    st.markdown("---")
    st.subheader("🧬 B-Cell Specific Analysis (S1B Signatures)")
    
    b_cell_types = [c for c in plot_df[CT_COL].dropna().unique() if 'b cell' in str(c).lower() or 'b' in str(c).lower().split() or 'plasma' in str(c).lower()]
    
    if b_cell_types:
        b_df = plot_df[plot_df[CT_COL].isin(b_cell_types)].copy()
        tb1, tb2, tb3, tb4, tb5 = st.tabs(["Composition", "Viral Load", "Leiden Clusters", "Subtype Markers (All Cells)", "Subtype Markers (UMAP)"])
        
        with tb1:
            st.plotly_chart(px.histogram(b_df, x=CT_COL, color='Disease_Condition (Detail)', barnorm='percent', title="B-Cell Subtypes by Disease Condition"), use_container_width=True)
            
        with tb2:
            b_infected = b_df[b_df['Viral_Counts'] >= 2]
            if not b_infected.empty: st.plotly_chart(px.box(b_infected, x=CT_COL, y='Viral_Counts', color='Disease_Condition (Detail)', points='all', title="Viral Load in Infected B-Cells"), use_container_width=True)
                
        with tb3:
            if 'leiden' in b_df.columns and 'Disease_Condition (Detail)' in b_df.columns:
                b_clus_dis = b_df.groupby(['leiden', 'Disease_Condition (Detail)']).size().reset_index(name='count')
                st.plotly_chart(px.bar(b_clus_dis, x='leiden', y='count', color='Disease_Condition (Detail)', title="Total B-Cell Counts per Cluster", barmode='stack'), use_container_width=True)
                
        with tb4:
            st.write("### 1. Global View: Static DotPlot (All Cell Types)")
            if BCELL_MARKERS:
                filtered_bcell_markers = {c: [g for g in genes if g in adata.var_names] for c, genes in BCELL_MARKERS.items() if [g for g in genes if g in adata.var_names]}
                
                if filtered_bcell_markers:
                    with st.spinner("Generating static DotPlot..."):
                        all_cell_types = sorted(adata.obs[CT_COL].dropna().unique())
                        
                        plot_height = min(15, max(5, len(all_cell_types) * 0.4))
                        plot_width = max(12, len(filtered_bcell_markers) * 2)
                        
                        fig, ax = plt.subplots(figsize=(plot_width, plot_height))
                        sc.pl.dotplot(adata, filtered_bcell_markers, groupby=CT_COL, dendrogram=False, ax=ax, show=False)
                        st.pyplot(fig)
                        
                    gene_to_cluster = {}
                    for cluster, genes in filtered_bcell_markers.items():
                        for g in genes:
                            if g not in gene_to_cluster:
                                gene_to_cluster[g] = cluster
                    
                    subtypes = list(filtered_bcell_markers.keys())
                    palette = px.colors.qualitative.Plotly + px.colors.qualitative.D3
                    cluster_color_map = {st_name: palette[i % len(palette)] for i, st_name in enumerate(subtypes)}

                    def format_gene_label(g):
                        cluster = gene_to_cluster.get(g, '')
                        color = cluster_color_map.get(cluster, 'black')
                        return f"<span style='color:{color}'>{g} ({cluster})</span>"
                    
                    flat_genes = []
                    for genes in filtered_bcell_markers.values():
                        for g in genes:
                            if g not in flat_genes: flat_genes.append(g)
                    
                    ordered_gene_labels = [format_gene_label(x) for x in flat_genes]
                            
                    with st.expander("📥 View & Download General DotPlot Data"):
                        df_expr = sc.get.obs_df(adata, keys=flat_genes + [CT_COL])
                        mean_df = df_expr.groupby(CT_COL).mean().reset_index()
                        
                        gene_cols = [c for c in df_expr.columns if c != CT_COL]
                        frac_df = df_expr.groupby(CT_COL)[gene_cols].apply(lambda x: (x > 0).mean()).reset_index()
                        
                        mean_melt = mean_df.melt(id_vars=CT_COL, var_name='Gene', value_name='Mean_Expression')
                        frac_melt = frac_df.melt(id_vars=CT_COL, var_name='Gene', value_name='Fraction_Expressing')
                        
                        export_df = pd.merge(mean_melt, frac_melt, on=[CT_COL, 'Gene'])
                        export_df['Gene'] = export_df['Gene'].apply(lambda x: f"{x} ({gene_to_cluster.get(x, '')})")
                        
                        st.dataframe(export_df)
                        
                        csv_data = export_df.to_csv(index=False).encode('utf-8')
                        st.download_button(
                            label="Download General Data (CSV)",
                            data=csv_data,
                            file_name="bcell_markers_expression.csv",
                            mime="text/csv",
                            key="dl_general"
                        )
                        
                    st.markdown("---")
                    st.write("### 2. Interactive Subtype Comparisons")
                    
                    cell_type_sel = st.multiselect("Select Cell Types for comparison:", all_cell_types, default=all_cell_types)
                    
                    if not cell_type_sel:
                        st.warning("⚠️ Please select at least one cell type from the dropdown to view comparisons.")
                    else:
                        required_obs = [CT_COL]
                        if 'Disease_Condition (Detail)' in adata.obs.columns: required_obs.append('Disease_Condition (Detail)')
                        if 'Infection' in adata.obs.columns: required_obs.append('Infection')
                        
                        df_expr_all = sc.get.obs_df(adata, keys=flat_genes + required_obs)
                        
                        # A. Disease Condition Comparison
                        if 'Disease_Condition (Detail)' in required_obs:
                            st.markdown("<br>", unsafe_allow_html=True)
                            st.write("#### Disease Condition (Healthy vs MS)")
                            
                            dis_mode = st.radio("Disease Comparison Mode:", ["Detailed (HS, Active, Stable)", "Combined (HS vs MS)"], horizontal=True, key="dis_mode")
                            
                            group_col = 'Disease_Condition (Detail)' if 'Detailed' in dis_mode else 'Disease_Group'
                            
                            if group_col == 'Disease_Group':
                                df_expr_all['Disease_Group'] = df_expr_all['Disease_Condition (Detail)'].apply(
                                    lambda x: 'MS' if 'Active' in str(x) or 'Stable' in str(x) else ('HS' if 'HS' in str(x) or 'Healthy' in str(x) else str(x))
                                )
                            
                            cols_to_drop = [c for c in ['Infection', 'Disease_Condition (Detail)', 'Disease_Group'] if c in df_expr_all.columns and c != group_col]
                            grouped_cond = df_expr_all.drop(columns=cols_to_drop).groupby([CT_COL, group_col])
                            
                            mean_df_cond = grouped_cond.mean().reset_index()
                            frac_df_cond = grouped_cond[gene_cols].apply(lambda x: (x > 0).mean()).reset_index()
                            
                            mean_cond_melt = mean_df_cond.melt(id_vars=[CT_COL, group_col], var_name='Gene', value_name='Mean_Expression')
                            frac_cond_melt = frac_df_cond.melt(id_vars=[CT_COL, group_col], var_name='Gene', value_name='Fraction_Expression')
                            merged_cond = pd.merge(mean_cond_melt, frac_cond_melt, on=[CT_COL, group_col, 'Gene'])
                            
                            df_plot_cond = merged_cond[merged_cond[CT_COL].isin(cell_type_sel)].copy()
                            
                            if not df_plot_cond.empty:
                                def color_cond(c):
                                    s = str(c)
                                    if 'HS' in s or 'Healthy' in s: return f"<b><span style='color:#2CA02C'>{s}</span></b>"
                                    if 'Active' in s: return f"<b><span style='color:#D62728'>{s}</span></b>"
                                    if 'Stable' in s: return f"<b><span style='color:#1F77B4'>{s}</span></b>"
                                    if s == 'MS': return f"<b><span style='color:#FF7F0E'>{s}</span></b>"
                                    return s
    
                                df_plot_cond['Condition_Colored'] = df_plot_cond[group_col].apply(color_cond)
                                df_plot_cond['Gene_Annotated'] = df_plot_cond['Gene'].apply(format_gene_label)
                                df_plot_cond['Y_Composite'] = df_plot_cond[CT_COL].astype(str) + " - " + df_plot_cond['Condition_Colored'].astype(str)
                                
                                plot_height_cond = max(400, len(cell_type_sel) * 80)
    
                                fig_cond = px.scatter(
                                    df_plot_cond, x='Gene_Annotated', y='Y_Composite', size='Fraction_Expression',
                                    color='Mean_Expression', color_continuous_scale='Reds', height=plot_height_cond, size_max=10, 
                                    hover_data={CT_COL: True, 'Y_Composite': False, group_col: True, 'Gene_Annotated': False, 'Gene': True}
                                )
                                
                                fig_cond.update_xaxes(categoryorder='array', categoryarray=ordered_gene_labels, tickmode='linear', tickfont=dict(size=10))
                                fig_cond.update_yaxes(categoryorder='category descending', title="")
                                fig_cond.update_layout(xaxis_title="", xaxis_tickangle=-90, margin=dict(l=0, r=0, t=10, b=180))
                                st.plotly_chart(fig_cond, use_container_width=True)
                                
                                with st.expander("📥 View & Download Disease Comparison Data"):
                                    pivot_cond = merged_cond.pivot_table(index=[CT_COL, 'Gene'], columns=group_col, values=['Mean_Expression', 'Fraction_Expression'])
                                    pivot_cond.columns = [f"{col[0]}_{col[1]}" for col in pivot_cond.columns]
                                    pivot_cond = pivot_cond.reset_index()
                                    pivot_cond['Gene'] = pivot_cond['Gene'].apply(lambda x: f"{x} ({gene_to_cluster.get(x, '')})")
                                    
                                    display_cond_table = pivot_cond[pivot_cond[CT_COL].isin(cell_type_sel)]
                                    st.dataframe(display_cond_table)
                                    csv_data_cond = display_cond_table.to_csv(index=False).encode('utf-8')
                                    st.download_button(label="Download Disease Comparison Table (CSV)", data=csv_data_cond, file_name="bcell_markers_condition_pivot.csv", mime="text/csv", key="dl_cond_pivot")
                            else:
                                st.warning("No condition data available for the selected cell types.")
    
                        # B. Infection Comparison
                        if 'Infection' in required_obs:
                            st.markdown("<br>", unsafe_allow_html=True)
                            st.write("#### Infection Comparison")
                            
                            if 'Disease_Condition (Detail)' in df_expr_all.columns and 'Disease_Group' not in df_expr_all.columns:
                                df_expr_all['Disease_Group'] = df_expr_all['Disease_Condition (Detail)'].apply(
                                    lambda x: 'MS' if 'Active' in str(x) or 'Stable' in str(x) else ('HC' if 'HC' in str(x) or 'Healthy' in str(x) else str(x))
                                )
                                
                            inf_mode_opts = ["Mock vs EBV (All)"]
                            if 'Disease_Group' in df_expr_all.columns:
                                inf_mode_opts.extend(["HC Only (Mock vs EBV)", "MS Only (Mock vs EBV)"])
                                
                            inf_mode = st.radio("Infection Comparison Mode:", inf_mode_opts, horizontal=True, key="inf_mode")
                            
                            if inf_mode == "Mock vs EBV (All)":
                                cols_to_drop = [c for c in ['Disease_Condition (Detail)', 'Disease_Group'] if c in df_expr_all.columns]
                                sub_df = df_expr_all.drop(columns=cols_to_drop)
                            elif "HC Only" in inf_mode:
                                sub_df = df_expr_all[df_expr_all['Disease_Group'] == 'HC'].copy()
                                cols_to_drop = [c for c in ['Disease_Condition (Detail)', 'Disease_Group'] if c in sub_df.columns]
                                sub_df = sub_df.drop(columns=cols_to_drop)
                            elif "MS Only" in inf_mode:
                                sub_df = df_expr_all[df_expr_all['Disease_Group'] == 'MS'].copy()
                                cols_to_drop = [c for c in ['Disease_Condition (Detail)', 'Disease_Group'] if c in sub_df.columns]
                                sub_df = sub_df.drop(columns=cols_to_drop)
                                
                            grouped_inf = sub_df.groupby([CT_COL, 'Infection'])
                            y_col = 'Infection'
                                
                            mean_df_inf = grouped_inf.mean().reset_index()
                            frac_df_inf = grouped_inf[gene_cols].apply(lambda x: (x > 0).mean()).reset_index()
                            
                            mean_inf_melt = mean_df_inf.melt(id_vars=[CT_COL, y_col], var_name='Gene', value_name='Mean_Expression')
                            frac_inf_melt = frac_df_inf.melt(id_vars=[CT_COL, y_col], var_name='Gene', value_name='Fraction_Expression')
                            merged_inf = pd.merge(mean_inf_melt, frac_inf_melt, on=[CT_COL, y_col, 'Gene'])
                            
                            pivot_inf = merged_inf.pivot_table(index=[CT_COL, 'Gene'], columns=y_col, values=['Mean_Expression', 'Fraction_Expression']).fillna(0)
                            pivot_inf.columns = [f"{col[0]}_{col[1]}" for col in pivot_inf.columns]
                            pivot_inf = pivot_inf.reset_index()
                            
                            df_plot_inf = merged_inf[merged_inf[CT_COL].isin(cell_type_sel)].copy()
                            
                            if not df_plot_inf.empty:
                                def color_inf(c):
                                    s = str(c)
                                    if 'Mock' in s: return f"<b><span style='color:#7F7F7F'>{s}</span></b>"
                                    return f"<b><span style='color:#9467BD'>{s}</span></b>"
    
                                df_plot_inf['Y_Colored'] = df_plot_inf[y_col].apply(color_inf)
                                df_plot_inf['Gene_Annotated'] = df_plot_inf['Gene'].apply(format_gene_label)
                                df_plot_inf['Y_Composite'] = df_plot_inf[CT_COL].astype(str) + " - " + df_plot_inf['Y_Colored'].astype(str)
                                
                                plot_height_inf = max(400, len(cell_type_sel) * 80)
    
                                fig_inf = px.scatter(
                                    df_plot_inf, x='Gene_Annotated', y='Y_Composite', size='Fraction_Expression',
                                    color='Mean_Expression', color_continuous_scale='Reds', height=plot_height_inf, size_max=10, 
                                    hover_data={CT_COL: True, 'Y_Composite': False, y_col: True, 'Gene_Annotated': False, 'Gene': True}
                                )
                                
                                fig_inf.update_xaxes(categoryorder='array', categoryarray=ordered_gene_labels, tickmode='linear', tickfont=dict(size=10))
                                fig_inf.update_yaxes(categoryorder='category descending', title="")
                                fig_inf.update_layout(xaxis_title="", xaxis_tickangle=-90, margin=dict(l=0, r=0, t=10, b=180))
                                st.plotly_chart(fig_inf, use_container_width=True)
                                
                                with st.expander("📥 View & Download Infection Comparison Data"):
                                    display_table = pivot_inf[pivot_inf[CT_COL].isin(cell_type_sel)].copy()
                                    display_table['Gene'] = display_table['Gene'].apply(lambda x: f"{x} ({gene_to_cluster.get(x, '')})")
                                    
                                    st.dataframe(display_table)
                                    csv_data_inf = display_table.to_csv(index=False).encode('utf-8')
                                    st.download_button(label="Download Infection Comparison Table (CSV)", data=csv_data_inf, file_name="bcell_multi_infection_pivot.csv", mime="text/csv", key="dl_inf_pivot")
                            else:
                                st.warning("No valid data available for the selected cell types under this comparison mode.")
                            
                else:
                    st.warning("No B-cell markers from the CSV were found in adata.var_names.")
            else:
                st.warning("B-cell markers CSV not loaded.")

        with tb5:
            st.write("Static UMAP Expression of Top B-Cell Subtype Markers")
            if BCELL_MARKERS:
                all_b_genes = []
                for cluster, genes in BCELL_MARKERS.items():
                    all_b_genes.extend(genes)
                
                valid_b_genes = sorted(list(set([g for g in all_b_genes if g in adata.var_names])))
                
                if valid_b_genes:
                    selected_umap_genes = st.multiselect("Select genes to overlay on UMAP:", valid_b_genes, default=valid_b_genes[:4])
                    
                    if selected_umap_genes:
                        with st.spinner("Generating static UMAPs. This may take a moment..."):
                            fig = sc.pl.umap(adata, color=selected_umap_genes, cmap='Reds', show=False, return_fig=True)
                            if fig is None: fig = plt.gcf()
                            st.pyplot(fig)
                            plt.close(fig)
                else:
                    st.warning("No valid B-cell markers found in dataset.")
            else:
                st.warning("B-cell markers CSV not loaded.")

if color_by == 'Disease Comparison':
    st.markdown("---")
    st.subheader("🏥 MS Markers across Disease Conditions")
    valid_ms = [g for g in MS_MARKER_LIST if g in adata.var_names]
    if valid_ms and 'Disease_Condition (Detail)' in adata.obs.columns:
        conditions = sorted([str(c) for c in plot_df['Disease_Condition (Detail)'].dropna().unique() if str(c) != 'nan'])
        if conditions:
            sub_adata = adata[plot_df.index]
            umap_coords = sub_adata.obsm['X_umap']
            with st.spinner("Generating static UMAPs for all MS Markers. This might take a moment..."):
                for marker in valid_ms:
                    fig, axes = plt.subplots(1, len(conditions), figsize=(4 * len(conditions), 4))
                    if len(conditions) == 1: axes = [axes]
                    expr = sub_adata[:, marker].X
                    expr_flat = expr.toarray().flatten() if sp.issparse(expr) else np.asarray(expr).flatten()
                    vmax = np.max(expr_flat) if np.max(expr_flat) > 0 else 1.0 
                    
                    for i, cond in enumerate(conditions):
                        ax = axes[i]
                        ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c='#e0e0e0', s=1, alpha=0.5, edgecolors='none')
                        cond_mask = (sub_adata.obs['Disease_Condition (Detail)'] == cond).values
                        cond_coords = umap_coords[cond_mask]
                        cond_expr = expr_flat[cond_mask]
                        sort_idx = np.argsort(cond_expr)
                        ax.scatter(cond_coords[sort_idx, 0], cond_coords[sort_idx, 1], c=cond_expr[sort_idx], cmap='Reds', s=3, vmin=0, vmax=vmax, edgecolors='none')
                        ax.set_title(f"{marker} ({cond})", fontsize=12)
                        ax.axis('off')
                        
                    plt.tight_layout()
                    st.pyplot(fig)
                    plt.close(fig)

if color_by == 'Mock Breakdown':
    st.markdown("---")
    st.subheader("🧪 Mock Infection Breakdown")
    if 'Infection' in plot_df.columns and 'Disease_Condition (Detail)' in plot_df.columns and 'Day' in plot_df.columns:
        mock_df = plot_df[plot_df['Infection'].astype(str).str.contains('Mock', case=False, na=False)].copy()
        if mock_df.empty:
            st.warning("No cells labeled as 'Mock' found in the 'Infection' column.")
        else:
            mock_df['Day_Str'] = mock_df['Day'].apply(lambda x: f"Day {int(get_day_sort_key(x))}" if get_day_sort_key(x) != 9999 else str(x))
            cond_strs = mock_df['Disease_Condition (Detail)'].astype(str).tolist()
            day_strs = mock_df['Day_Str'].tolist()
            mock_df['Disease_Day'] = [f"{c} - {d}" for c, d in zip(cond_strs, day_strs)]
            
            mock_df['Day_Sorted'] = mock_df['Day'].astype(str).apply(get_day_sort_key).astype(int)
            mock_df = mock_df.sort_values(by=['Disease_Condition (Detail)', 'Day_Sorted'])
            
            fig = px.histogram(mock_df, x='Disease_Day', color=CT_COL, barnorm='percent', title="Mock Composition by Disease & Day")
            fig.update_xaxes(categoryorder='array', categoryarray=mock_df['Disease_Day'].unique())
            st.plotly_chart(fig, use_container_width=True)

if color_by == 'Trajectory Analysis':
    st.markdown("---")
    st.subheader("🛤️ Time-Aware Trajectory Analysis")
    
    if 'Day' not in adata.obs.columns:
        st.error("Column 'Day' not found in metadata. Cannot perform temporal trajectory analysis.")
    else:
        pseudo_col = 'dpt_pseudotime' if 'dpt_pseudotime' in adata.obs.columns else 'Pseudotime' if 'Pseudotime' in adata.obs.columns else None
        
        if pseudo_col:
            tt1, tt2, tt3, tt4 = st.tabs(["Temporal PAGA Graph", "Diffusion Map", "Pseudotime vs Actual Day", "Cell Types Across Time"])
            
            with tt1:
                st.write("**PAGA Connectivity Map (Clusters arranged chronologically by Pseudotime)**")
                
                try:
                    med_dpt = adata.obs.groupby('leiden')[pseudo_col].median()
                    med_dpt = med_dpt.reindex(adata.obs['leiden'].cat.categories).fillna(0).values
                    
                    sc.pl.paga(adata, show=False)
                    pos = adata.uns['paga']['pos'].copy()
                    
                    dpt_min, dpt_max = med_dpt.min(), med_dpt.max()
                    dpt_norm = (med_dpt - dpt_min) / (dpt_max - dpt_min + 1e-9)
                    
                    y_range = pos[:, 1].max() - pos[:, 1].min()
                    if y_range == 0: y_range = 1.0
                    
                    pos[:, 0] = dpt_norm * (y_range * 2.0)
                    
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
                    
                    sc.pl.paga(adata, pos=pos, color='leiden', ax=ax1, show=False, title="Clusters Ordered by Pseudotime")
                    sc.pl.paga(adata, pos=pos, color=pseudo_col, cmap='viridis', ax=ax2, show=False, title="Pseudotime Flow Gradient")
                    
                    plt.tight_layout()
                    st.pyplot(fig)
                except Exception as e:
                    st.warning(f"Could not generate ordered PAGA layout: {e}")
                    
            with tt2:
                st.write("**Diffusion Map (DC1 vs DC2) mapping Actual Time onto Trajectory Components**")
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
                sc.pl.diffmap(adata, color='Day', ax=ax1, show=False)
                sc.pl.diffmap(adata, color=pseudo_col, ax=ax2, show=False, cmap='viridis')
                plt.tight_layout()
                st.pyplot(fig)
                
            with tt3:
                st.write("**Correlation between Computed Pseudotime and Actual Day**")
                plot_df['Day_Sorted'] = plot_df['Day'].astype(str).apply(get_day_sort_key).astype(int)
                temp_df = plot_df.sort_values(by='Day_Sorted')
                fig_day = px.box(temp_df, x='Day', y=pseudo_col, color='Day', points='all', title=f"Pseudotime Distribution over Actual Days")
                st.plotly_chart(fig_day, use_container_width=True)
                
            with tt4:
                st.write("**Cellular Progression Profile**")
                fig_box = px.box(plot_df, x=CT_COL, y=pseudo_col, color='Day', title="Pseudotime Distribution by Cell Type (Split by Day)")
                fig_box.update_layout(xaxis={'categoryorder': 'median ascending'})
                st.plotly_chart(fig_box, use_container_width=True)
        else:
            st.info("Pseudotime could not be calculated automatically. Check dataset contents.")

@st.cache_data
def translate_markers(gene_list):
    mg = mygene.MyGeneInfo()
    clean_ids = [g.split('.')[0] if isinstance(g, str) and g.startswith('ENSG') else g for g in gene_list]
    ensg_only = list(set([g for g in clean_ids if str(g).startswith('ENSG')]))
    
    mapping = {orig: orig for orig in gene_list}
    if not ensg_only: return mapping
    
    try:
        res = mg.querymany(ensg_only, scopes='ensembl.gene', fields='symbol', species='human', as_dataframe=True, verbose=False)
        for orig, clean in zip(gene_list, clean_ids):
            if clean in res.index and 'symbol' in res.columns:
                sym = res.loc[clean, 'symbol']
                if isinstance(sym, pd.Series): sym = sym.iloc[0]
                if pd.notna(sym):
                    mapping[orig] = str(sym).upper()
        return mapping
    except:
        return mapping

if color_by in ['leiden', CT_COL]:
    st.markdown("---")
    st.subheader(f"🧬 Diagnostics: {color_by}")
    
    t_clus1, t_clus2 = st.tabs(["Disease Proportions", "Top & Bottom Markers"])
    
    with t_clus1:
        if 'Disease_Condition (Detail)' in plot_df.columns:
            df_clus_dis = plot_df.groupby([color_by, 'Disease_Condition (Detail)']).size().reset_index(name='count')
            df_clus_dis['percentage'] = df_clus_dis.groupby(color_by)['count'].transform(lambda x: x / x.sum() * 100)
            fig_bar = px.bar(df_clus_dis, x=color_by, y='percentage', color='Disease_Condition (Detail)', title=f"Disease Proportion per {color_by} (%)", barmode='stack')
            st.plotly_chart(fig_bar, use_container_width=True)
        else:
            st.info("'Disease_Condition (Detail)' not found.")

    with t_clus2:
        marker_key = f'markers_{color_by}'
        if marker_key not in adata.uns and 'rank_genes_groups' in adata.uns:
            marker_key = 'rank_genes_groups'

        if marker_key in adata.uns:
            try:
                # ---> BUG FIX: Pull the original cell type names directly from the marker dictionary <---
                # This ensures we capture pre-relabeled types like 'ILC3' and 'Tem/Trm cytotoxic T cells'
                original_groups = list(adata.uns[marker_key]['names'].dtype.names)
                original_groups = sorted([str(g) for g in original_groups if str(g) != 'nan'])
                
                all_marker_genes = []
                for grp in original_groups:
                    all_marker_genes.extend(adata.uns[marker_key]['names'][grp][:10])
                    all_marker_genes.extend(adata.uns[marker_key]['names'][grp][-10:])
                
                with st.spinner("Translating Ensembl IDs to Gene Symbols..."):
                    trans_dict = translate_markers(list(set(all_marker_genes)))
                
                cols = st.columns(3)
                for i, grp in enumerate(original_groups):
                    up = adata.uns[marker_key]['names'][grp][:10]
                    up_s = adata.uns[marker_key]['scores'][grp][:10]
                    down = adata.uns[marker_key]['names'][grp][-10:]
                    down_s = adata.uns[marker_key]['scores'][grp][-10:]
                    
                    up_mapped = [trans_dict.get(g, g) for g in up]
                    down_mapped = [trans_dict.get(g, g) for g in down]
                    
                    df_m = pd.DataFrame({'Gene': list(up_mapped) + list(down_mapped), 'Score': list(up_s) + list(down_s)}).sort_values('Score')
                    df_m['Original_ENSG'] = list(up) + list(down)
                    
                    # Highlight if it's one of the original complex names
                    title_str = f"{grp}"
                    
                    fig = px.bar(df_m, x='Score', y='Gene', orientation='h', title=title_str, color='Score', color_continuous_scale='RdBu_r', height=600, hover_data=['Original_ENSG'])
                    cols[i%3].plotly_chart(fig, use_container_width=True)
            except Exception as e: 
                st.error(f"Error plotting markers: {e}")
        else:
            st.warning(f"Marker genes have not been calculated for '{color_by}'. Run the add_celltype_markers.py script.")

if color_by == 'QC Metrics':
    st.markdown("---")
    if sensitivity_df is not None: st.table(sensitivity_df)
