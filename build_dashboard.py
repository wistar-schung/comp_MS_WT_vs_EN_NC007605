import os

html_content = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EBV Multi-Compartment Dashboard</title>
    <style>
        body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; display: flex; height: 100vh; overflow: hidden; background-color: #f8f9fa; }
        .sidebar { width: 400px; background-color: #ffffff; border-right: 1px solid #dee2e6; display: flex; flex-direction: column; overflow-y: auto; box-shadow: 2px 0 5px rgba(0,0,0,0.05); z-index: 10; flex-shrink: 0; }
        .header { padding: 25px 20px; border-bottom: 1px solid #dee2e6; background-color: #f8f9fa; position: sticky; top: 0; z-index: 20; }
        .header h2 { margin: 0; font-size: 1.3rem; color: #212529; }
        .header p { margin: 5px 0 0 0; font-size: 0.85rem; color: #6c757d; }
        .nav-section { padding: 20px 20px 8px 20px; font-weight: 700; font-size: 0.75rem; color: #adb5bd; text-transform: uppercase; letter-spacing: 1px; border-top: 1px solid #f1f3f5; margin-top: 10px; }
        .nav-section:first-of-type { border-top: none; margin-top: 0; }
        .nav-item { padding: 8px 20px 8px 30px; cursor: pointer; color: #495057; text-decoration: none; display: block; border-left: 3px solid transparent; font-size: 0.85rem; transition: all 0.2s ease; }
        .nav-item:hover { background-color: #f1f3f5; color: #212529; }
        .nav-item.active { background-color: #e7f1ff; color: #0d6efd; border-left-color: #0d6efd; font-weight: 600; }
        .content { flex-grow: 1; display: flex; flex-direction: column; background-color: #ffffff; overflow: hidden; }
        iframe { flex-grow: 1; border: none; width: 100%; height: 100%; }
        .sidebar::-webkit-scrollbar { width: 6px; }
        .sidebar::-webkit-scrollbar-track { background: #f1f1f1; }
        .sidebar::-webkit-scrollbar-thumb { background: #c1c1c1; border-radius: 3px; }
    </style>
</head>
<body>
    <div class="sidebar">
        <div class="header">
            <h2>🦠 EBV Interactomics</h2>
            <p>Master Analytics Dashboard</p>
        </div>
        
        <div class="nav-section">🌍 1. Global Overviews</div>
        <a class="nav-item active" onclick="loadPlot('Dashboards_01_Global/01_Dataset.html', this)">Dataset Breakdown (WT vs Enr)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_01_Global/02_CellTypes.html', this)">Global Cell Types</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_01_Global/03a_Global_Composition_Bar_Abs.html', this)">Global Composition Bar (Absolute)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_01_Global/03b_Global_Composition_Bar_Pct.html', this)">Global Composition Bar (Percentage)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_01_Global/04_Global_ViralLoad_Log.html', this)">Global Viral Load (Log1p)</a>
        
        <div class="nav-section">🧫 2. Compartment Overviews</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Bcell_00_Subtypes.html', this)">B-Cells: Subtypes Map</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Bcell_01_TargetScore.html', this)">B-Cells: ABC Signature Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Bcell_02_IFN_Score.html', this)">B-Cells: IFN Response Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Bcell_03_Timecourse_Score.html', this)">B-Cells: Timecourse Dynamics</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Tcell_00_Subtypes.html', this)">T-Cells: Subtypes Map</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Tcell_01_TargetScore.html', this)">T-Cells: Exhaustion Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Tcell_02_IFN_Score.html', this)">T-Cells: IFN Response Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Tcell_03_Timecourse_Score.html', this)">T-Cells: Timecourse Dynamics</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Myeloid_00_Subtypes.html', this)">Myeloid: Subtypes Map</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Myeloid_01_TargetScore.html', this)">Myeloid: Inflammatory Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Myeloid_02_IFN_Score.html', this)">Myeloid: IFN Response Score</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_02_Compartments/Myeloid_03_Timecourse_Score.html', this)">Myeloid: Timecourse Dynamics</a>

        <div class="nav-section">🧬 3. Pathway Dynamics (WT)</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_05_Pathway_Dynamics/Interactive_Pathway_Dynamics_WT.html', this)">Functional Enrichment Summary</a>

        <div class="nav-section">🎯 4. Target Cell Focus</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_06_CD16pos_NK_Focus/01_CD16pos_NK_Population_Dynamics.html', this)">CD16+ NK: WT Population Dynamics</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_06_CD16pos_NK_Focus/02_CD16pos_NK_EBV_Gene_Expression.html', this)">CD16+ NK: Enriched EBV Profile</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_06_ABC_Focus/01_ABC_Population_Dynamics.html', this)">ABC: WT Population Dynamics</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_06_ABC_Focus/02_ABC_EBV_Gene_Expression.html', this)">ABC: Enriched EBV Profile</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_06_NonClassicalMono_Focus/01_NonClassicalMono_Population_Dynamics.html', this)">Non-Class Mono: WT Population Dynamics</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_06_NonClassicalMono_Focus/02_NonClassicalMono_EBV_Gene_Expression.html', this)">Non-Class Mono: Enriched EBV Profile</a>

        <a class="nav-item" onclick="loadPlot('Dashboards_06_Naive_Helper_T_Focus/01_Naive_Helper_T_Population_Dynamics.html', this)">Naive Helper T: WT Population Dynamics</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_06_Naive_Helper_T_Focus/02_Naive_Helper_T_EBV_Gene_Expression.html', this)">Naive Helper T: Enriched EBV Profile</a>

        <div class="nav-section">📊 5. Differential Expression (WT)</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/01_CD16pos_NK_MS_vs_HC_Mock_D1_Volcano.html', this)">CD16+ NK: Active MS vs HC (D1 Mock)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/02_CD16pos_NK_ActiveMS_EBV_D7_vs_D1_Volcano.html', this)">CD16+ NK: Active MS EBV (D7 vs D1)</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/01_ABC_MS_vs_HC_Mock_D1_Volcano.html', this)">ABC: Active MS vs HC (D1 Mock)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/02_ABC_ActiveMS_EBV_D7_vs_D1_Volcano.html', this)">ABC: Active MS EBV (D7 vs D1)</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/01_NonClassicalMono_MS_vs_HC_Mock_D1_Volcano.html', this)">Non-Class Mono: Active MS vs HC (D1)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/02_NonClassicalMono_ActiveMS_EBV_D7_vs_D1_Volcano.html', this)">Non-Class Mono: Active MS EBV (D7 vs D1)</a>

        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/01_Naive_Helper_T_MS_vs_HC_Mock_D1_Volcano.html', this)">Naive Helper T: Active MS vs HC (D1)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_07_Differential_Expression/02_Naive_Helper_T_ActiveMS_EBV_D7_vs_D1_Volcano.html', this)">Naive Helper T: Active MS EBV (D7 vs D1)</a>

        <div class="nav-section">🔍 6. EBV Transcript Deep Dive</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Bcell_01_Top5_Genes_by_Disease.html', this)">B-Cells: Top 5 EBV Genes</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Bcell_02_Correlations_by_Disease.html', this)">B-Cells: EBV Drivers of ABC Score</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Tcell_01_Top5_Genes_by_Disease.html', this)">T-Cells: Top 5 EBV Genes</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Tcell_02_Correlations_by_Disease.html', this)">T-Cells: EBV Drivers of Exhaustion</a>
        
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Myeloid_01_Top5_Genes_by_Disease.html', this)">Myeloid: Top 5 EBV Genes</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_03_EBV_DeepDive/Myeloid_02_Correlations_by_Disease.html', this)">Myeloid: EBV Drivers of Inflammation</a>

        <div class="nav-section">📅 7. Day 1 Contrasts</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day1_UMAP_CellTypes.html', this)">Q0: UMAP Cell Types</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day1_UMAP_Disease.html', this)">Q0: UMAP Disease Groups</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day1_UMAP_Infection.html', this)">Q0: UMAP Infection Status</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day1_UMAP_ViralLoad.html', this)">Q0: UMAP Viral Load</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day1_WT_Viral_Load.html', this)">Q1: WT Viral Load (HC vs MS)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day1_WT_Composition_Shifts.html', this)">Q1: WT Cell Type Shifts</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q2_Q3_Day1_Enriched_EBV_DotPlot.html', this)">Q2/3: Enriched EBV DotPlot</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q4_Day1_Global_CellType_Predictions.html', this)">Q4: Global Cell Type Predictions</a>

        <div class="nav-section">📅 8. Day 7 Contrasts</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day7_UMAP_CellTypes.html', this)">Q0: UMAP Cell Types</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day7_UMAP_Disease.html', this)">Q0: UMAP Disease Groups</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day7_UMAP_Infection.html', this)">Q0: UMAP Infection Status</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day7_UMAP_ViralLoad.html', this)">Q0: UMAP Viral Load</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day7_WT_Viral_Load.html', this)">Q1: WT Viral Load (HC vs MS)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day7_WT_Composition_Shifts.html', this)">Q1: WT Cell Type Shifts</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q2_Q3_Day7_Enriched_EBV_DotPlot.html', this)">Q2/3: Enriched EBV DotPlot</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q4_Day7_Global_CellType_Predictions.html', this)">Q4: Global Cell Type Predictions</a>

        <div class="nav-section">📅 9. Day 15 Contrasts</div>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day15_UMAP_CellTypes.html', this)">Q0: UMAP Cell Types</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day15_UMAP_Disease.html', this)">Q0: UMAP Disease Groups</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day15_UMAP_Infection.html', this)">Q0: UMAP Infection Status</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q0_Day15_UMAP_ViralLoad.html', this)">Q0: UMAP Viral Load</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day15_WT_Viral_Load.html', this)">Q1: WT Viral Load (HC vs MS)</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q1_Day15_WT_Composition_Shifts.html', this)">Q1: WT Cell Type Shifts</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q2_Q3_Day15_Enriched_EBV_DotPlot.html', this)">Q2/3: Enriched EBV DotPlot</a>
        <a class="nav-item" onclick="loadPlot('Dashboards_04_Day_Contrasts/Q4_Day15_Global_CellType_Predictions.html', this)">Q4: Global Cell Type Predictions</a>

        <div style="height: 60px;"></div>
    </div>
    
    <div class="content">
        <iframe id="plot-frame" src="Dashboards_01_Global/01_Dataset.html"></iframe>
    </div>

    <script>
        function loadPlot(filepath, element) {
            document.getElementById('plot-frame').src = filepath;
            var items = document.getElementsByClassName('nav-item');
            for(var i=0; i<items.length; i++) {
                items[i].classList.remove('active');
            }
            if(element) {
                element.classList.add('active');
            }
        }
    </script>
</body>
</html>
"""

output_file = "Dashboard_Index.html"
with open(output_file, "w", encoding='utf-8') as f: 
    f.write(html_content)
    
print(f"✅ Dashboard Sidebar updated successfully! Double-click '{output_file}' to view your complete multi-compartment analysis.")