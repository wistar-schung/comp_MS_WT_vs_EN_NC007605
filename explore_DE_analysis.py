"""
Interactive DE Analysis Explorer

Generates an interactive HTML dashboard to explore all differential expression analyses.
Links volcano plots to their corresponding functional enrichment results.

Outputs:
- Master_DE_Explorer.html: Interactive dashboard with searchable/filterable interface
- Organized by cell type and contrast type
"""

import os
import json
import pandas as pd
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# ==========================================
# CONFIGURATION
# ==========================================
DIR_DE = "Dashboards_07_Differential_Expression_Pseudobulk"
DIR_DE_ENRICHMENT = "Dashboards_07_Differential_Expression_Pseudobulk/Functional_Enrichment"
OUTPUT_FILE = "Master_DE_Explorer.html"

CELL_TYPES = {
    'CD16pos_NK': 'CD16+ NK cells',
    'ABC': 'Age-associated B cells',
    'NonClassMono': 'Non-classical monocytes',
    'Naive_Helper_T': 'Tcm/Naive helper T cells'
}

CONTRAST_TYPES = {
    'Base': 'Baseline Disease Signatures',
    'Viral': 'Direct Viral Reaction',
    'Kin': 'Viral Kinetics',
    'Int': 'Statistical Interactions'
}

# ==========================================
# HELPER FUNCTIONS
# ==========================================

def categorize_contrast(prefix):
    """Categorize a contrast by its prefix."""
    for key, label in CONTRAST_TYPES.items():
        if prefix.startswith(key):
            return key, label
    return 'Other', 'Other'


def extract_files_by_celltype():
    """Scan directories and organize files by cell type and contrast."""
    data = defaultdict(lambda: defaultdict(list))
    
    if not os.path.exists(DIR_DE):
        print(f"Warning: {DIR_DE} not found")
        return data
    
    # Scan volcano plots
    for filename in os.listdir(DIR_DE):
        if filename.endswith('_Volcano.html'):
            # Parse filename: {CellType}_{Prefix}_Volcano.html
            parts = filename.replace('_Volcano.html', '').split('_')
            cell_type = parts[0]
            prefix = '_'.join(parts[1:])
            
            contrast_key, contrast_label = categorize_contrast(prefix)
            
            data[cell_type][contrast_key].append({
                'prefix': prefix,
                'volcano_file': filename,
                'volcano_path': os.path.join(DIR_DE, filename),
                'sig_genes_file': f"{cell_type}_{prefix}_Sig_Genes.csv",
                'contrast_label': contrast_label,
            })
    
    return data


def find_enrichment_files(cell_type, prefix):
    """Find enrichment files for a given cell type and contrast."""
    if not os.path.exists(DIR_DE_ENRICHMENT):
        return {}
    
    enrichment_files = {
        'upregulated': None,
        'downregulated': None
    }
    
    for direction in enrichment_files.keys():
        html_file = f"{cell_type}_{prefix}_{direction}_Enrichment.html"
        csv_file = f"{cell_type}_{prefix}_{direction}_Enrichment.csv"
        
        html_path = os.path.join(DIR_DE_ENRICHMENT, html_file)
        csv_path = os.path.join(DIR_DE_ENRICHMENT, csv_file)
        
        if os.path.exists(html_path):
            enrichment_files[direction] = {
                'html': html_file,
                'csv': csv_file if os.path.exists(csv_path) else None
            }
    
    return enrichment_files


def load_sig_genes_summary(cell_type, prefix):
    """Load summary stats from significant genes CSV."""
    csv_path = os.path.join(DIR_DE, f"{cell_type}_{prefix}_Sig_Genes.csv")
    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path)
            n_up = len(df[df['Significance'].str.contains('Hyper-|True', regex=True, na=False)])
            n_down = len(df[df['Significance'].str.contains('Failed|False', regex=True, na=False)])
            return {
                'total': len(df),
                'upregulated': n_up,
                'downregulated': n_down
            }
        except Exception as e:
            print(f"Warning: Could not read {csv_path}: {e}")
    return {'total': 0, 'upregulated': 0, 'downregulated': 0}


# ==========================================
# HTML GENERATION
# ==========================================

def generate_html_explorer(data):
    """Generate interactive HTML explorer."""
    
    html_parts = []
    
    # Header
    html_parts.append(f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>DE Analysis Explorer</title>
        <style>
            * {{
                margin: 0;
                padding: 0;
                box-sizing: border-box;
            }}
            
            body {{
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                min-height: 100vh;
                padding: 20px;
                color: #333;
            }}
            
            .container {{
                max-width: 1400px;
                margin: 0 auto;
                background: white;
                border-radius: 12px;
                box-shadow: 0 20px 60px rgba(0,0,0,0.3);
                overflow: hidden;
            }}
            
            .header {{
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 40px;
                text-align: center;
            }}
            
            .header h1 {{
                font-size: 2.5em;
                margin-bottom: 10px;
            }}
            
            .header p {{
                font-size: 1.1em;
                opacity: 0.9;
            }}
            
            .controls {{
                padding: 20px 40px;
                background: #f8f9fa;
                border-bottom: 1px solid #dee2e6;
                display: flex;
                gap: 20px;
                flex-wrap: wrap;
                align-items: center;
            }}
            
            .search-box {{
                flex: 1;
                min-width: 250px;
            }}
            
            .search-box input {{
                width: 100%;
                padding: 10px 15px;
                border: 2px solid #dee2e6;
                border-radius: 6px;
                font-size: 1em;
                transition: border-color 0.2s;
            }}
            
            .search-box input:focus {{
                outline: none;
                border-color: #667eea;
                box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
            }}
            
            .filters {{
                display: flex;
                gap: 10px;
                flex-wrap: wrap;
            }}
            
            .filter-btn {{
                padding: 8px 16px;
                border: 2px solid #dee2e6;
                background: white;
                border-radius: 6px;
                cursor: pointer;
                font-weight: 500;
                transition: all 0.2s;
            }}
            
            .filter-btn:hover {{
                border-color: #667eea;
                color: #667eea;
            }}
            
            .filter-btn.active {{
                background: #667eea;
                color: white;
                border-color: #667eea;
            }}
            
            .content {{
                padding: 40px;
            }}
            
            .cell-type-section {{
                margin-bottom: 50px;
            }}
            
            .cell-type-title {{
                font-size: 1.8em;
                color: #667eea;
                margin-bottom: 20px;
                padding-bottom: 10px;
                border-bottom: 3px solid #667eea;
            }}
            
            .contrast-type-group {{
                margin-bottom: 30px;
                padding: 20px;
                background: #f8f9fa;
                border-radius: 8px;
                border-left: 4px solid #764ba2;
            }}
            
            .contrast-type-label {{
                font-size: 1.2em;
                font-weight: 600;
                color: #764ba2;
                margin-bottom: 15px;
            }}
            
            .contrast-card {{
                background: white;
                border: 1px solid #dee2e6;
                border-radius: 8px;
                padding: 20px;
                margin-bottom: 15px;
                transition: all 0.2s;
            }}
            
            .contrast-card:hover {{
                box-shadow: 0 8px 16px rgba(0,0,0,0.1);
                border-color: #667eea;
            }}
            
            .contrast-header {{
                display: flex;
                justify-content: space-between;
                align-items: start;
                margin-bottom: 15px;
                flex-wrap: wrap;
                gap: 10px;
            }}
            
            .contrast-title {{
                font-size: 1.1em;
                font-weight: 600;
                color: #333;
            }}
            
            .contrast-stats {{
                display: flex;
                gap: 20px;
                font-size: 0.9em;
            }}
            
            .stat-badge {{
                padding: 4px 12px;
                background: #e9ecef;
                border-radius: 20px;
                font-weight: 500;
            }}
            
            .stat-badge.up {{
                background: #d62728;
                color: white;
            }}
            
            .stat-badge.down {{
                background: #1f77b4;
                color: white;
            }}
            
            .contrast-links {{
                display: flex;
                gap: 15px;
                flex-wrap: wrap;
            }}
            
            .link-group {{
                display: flex;
                flex-direction: column;
                gap: 8px;
            }}
            
            .link-label {{
                font-size: 0.85em;
                font-weight: 600;
                color: #666;
                text-transform: uppercase;
            }}
            
            .btn {{
                display: inline-flex;
                align-items: center;
                padding: 10px 16px;
                border: none;
                border-radius: 6px;
                font-weight: 500;
                cursor: pointer;
                transition: all 0.2s;
                text-decoration: none;
                font-size: 0.95em;
                gap: 6px;
            }}
            
            .btn-primary {{
                background: #667eea;
                color: white;
            }}
            
            .btn-primary:hover {{
                background: #5568d3;
                transform: translateY(-2px);
                box-shadow: 0 4px 12px rgba(102, 126, 234, 0.4);
            }}
            
            .btn-secondary {{
                background: #6c757d;
                color: white;
            }}
            
            .btn-secondary:hover {{
                background: #5a6268;
                transform: translateY(-2px);
                box-shadow: 0 4px 12px rgba(108, 117, 125, 0.4);
            }}
            
            .btn-success {{
                background: #28a745;
                color: white;
            }}
            
            .btn-success:hover {{
                background: #218838;
                transform: translateY(-2px);
                box-shadow: 0 4px 12px rgba(40, 167, 69, 0.4);
            }}
            
            .btn-sm {{
                padding: 6px 12px;
                font-size: 0.85em;
            }}
            
            .no-results {{
                text-align: center;
                padding: 40px;
                color: #999;
            }}
            
            .footer {{
                padding: 20px 40px;
                background: #f8f9fa;
                text-align: center;
                color: #666;
                font-size: 0.9em;
                border-top: 1px solid #dee2e6;
            }}
            
            .icon {{
                font-size: 1.2em;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1>🔬 Differential Expression Analysis Explorer</h1>
                <p>Interactive dashboard linking volcano plots to functional enrichment analyses</p>
            </div>
            
            <div class="controls">
                <div class="search-box">
                    <input type="text" id="searchInput" placeholder="Search by contrast name (e.g., 'Active', 'EBV', 'Day 15')...">
                </div>
                <div class="filters" id="filterButtons"></div>
            </div>
            
            <div class="content" id="content">
    """)
    
    # Generate content for each cell type
    cell_type_count = 0
    for cell_type, contrast_data in sorted(data.items()):
        if cell_type not in CELL_TYPES:
            continue
        
        cell_type_count += 1
        cell_type_name = CELL_TYPES.get(cell_type, cell_type)
        
        html_parts.append(f"""
            <div class="cell-type-section" data-cell-type="{cell_type}">
                <div class="cell-type-title">{cell_type_name}</div>
        """)
        
        # Group by contrast type
        for contrast_key in sorted(CONTRAST_TYPES.keys()):
            if contrast_key not in contrast_data:
                continue
            
            contrast_label = CONTRAST_TYPES[contrast_key]
            html_parts.append(f"""
                <div class="contrast-type-group" data-contrast-type="{contrast_key}">
                    <div class="contrast-type-label">{contrast_label}</div>
            """)
            
            # Add each contrast
            for contrast_info in contrast_data[contrast_key]:
                prefix = contrast_info['prefix']
                volcano_file = contrast_info['volcano_file']
                contrast_label_full = contrast_info['contrast_label']
                
                # Load summary stats
                stats = load_sig_genes_summary(cell_type, prefix)
                
                # Find enrichment files
                enrichment = find_enrichment_files(cell_type, prefix)
                
                # Format title from prefix
                title = prefix.replace('_', ' ')
                
                html_parts.append(f"""
                    <div class="contrast-card" data-contrast="{prefix}">
                        <div class="contrast-header">
                            <div class="contrast-title">{title}</div>
                            <div class="contrast-stats">
                                <span class="stat-badge">Total: {stats['total']}</span>
                                <span class="stat-badge up">↑ {stats['upregulated']}</span>
                                <span class="stat-badge down">↓ {stats['downregulated']}</span>
                            </div>
                        </div>
                        <div class="contrast-links">
                """)
                
                # Volcano plot link
                html_parts.append(f"""
                    <div class="link-group">
                        <div class="link-label">📊 Volcano Plot</div>
                        <a href="{volcano_file}" class="btn btn-primary" target="_blank">
                            <span class="icon">📈</span>View Plot
                        </a>
                    </div>
                """)
                
                # Enrichment links
                if enrichment.get('upregulated'):
                    html_parts.append(f"""
                    <div class="link-group">
                        <div class="link-label">🔼 Upregulated</div>
                        <a href="{os.path.join('Functional_Enrichment', enrichment['upregulated']['html'])}" class="btn btn-success btn-sm" target="_blank">
                            Enrichment Plot
                        </a>
                    """)
                    if enrichment['upregulated']['csv']:
                        html_parts.append(f"""
                        <a href="{os.path.join('Functional_Enrichment', enrichment['upregulated']['csv'])}" class="btn btn-secondary btn-sm" target="_blank">
                            Data (CSV)
                        </a>
                    """)
                    html_parts.append("</div>")
                
                if enrichment.get('downregulated'):
                    html_parts.append(f"""
                    <div class="link-group">
                        <div class="link-label">🔽 Downregulated</div>
                        <a href="{os.path.join('Functional_Enrichment', enrichment['downregulated']['html'])}" class="btn btn-success btn-sm" target="_blank">
                            Enrichment Plot
                        </a>
                    """)
                    if enrichment['downregulated']['csv']:
                        html_parts.append(f"""
                        <a href="{os.path.join('Functional_Enrichment', enrichment['downregulated']['csv'])}" class="btn btn-secondary btn-sm" target="_blank">
                            Data (CSV)
                        </a>
                    """)
                    html_parts.append("</div>")
                
                html_parts.append("""
                        </div>
                    </div>
                """)
            
            html_parts.append("""
                </div>
            """)
        
        html_parts.append("""
            </div>
        """)
    
    # Footer
    html_parts.append(f"""
            </div>
            
            <div class="footer">
                <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p>Total: {cell_type_count} cell types | Filter and search above to explore analyses</p>
            </div>
        </div>
        
        <script>
            // Search and filter functionality
            const searchInput = document.getElementById('searchInput');
            const filterButtons = document.getElementById('filterButtons');
            const contentDiv = document.getElementById('content');
            
            // Create filter buttons
            const cellTypes = {CELL_TYPES};
            const contrastTypes = {CONTRAST_TYPES};
            
            // Add filter buttons for contrast types
            Object.entries(contrastTypes).forEach(([key, label]) => {{
                const btn = document.createElement('button');
                btn.className = 'filter-btn active';
                btn.textContent = label;
                btn.dataset.filter = key;
                btn.onclick = () => applyFilters();
                filterButtons.appendChild(btn);
            }});
            
            function applyFilters() {{
                const searchTerm = searchInput.value.toLowerCase();
                const activeFilters = Array.from(document.querySelectorAll('.filter-btn.active')).map(b => b.dataset.filter);
                
                // Show/hide cell type sections
                document.querySelectorAll('.cell-type-section').forEach(section => {{
                    let sectionVisible = false;
                    
                    // Show/hide contrast groups
                    section.querySelectorAll('.contrast-type-group').forEach(group => {{
                        const contrastType = group.dataset.contrastType;
                        const isTypeActive = activeFilters.includes(contrastType);
                        
                        let groupVisible = false;
                        
                        // Show/hide individual contrasts
                        group.querySelectorAll('.contrast-card').forEach(card => {{
                            const contrastText = card.dataset.contrast.toLowerCase();
                            const matches = searchTerm === '' || contrastText.includes(searchTerm);
                            
                            card.style.display = (matches && isTypeActive) ? 'block' : 'none';
                            if (matches && isTypeActive) {{
                                groupVisible = true;
                                sectionVisible = true;
                            }}
                        }});
                        
                        group.style.display = groupVisible ? 'block' : 'none';
                    }});
                    
                    section.style.display = sectionVisible ? 'block' : 'none';
                }});
            }}
            
            // Filter button toggle
            document.querySelectorAll('.filter-btn').forEach(btn => {{
                btn.addEventListener('click', (e) => {{
                    e.target.classList.toggle('active');
                    applyFilters();
                }});
            }});
            
            // Search input
            searchInput.addEventListener('input', applyFilters);
            
            // Initialize
            applyFilters();
        </script>
    </body>
    </html>
    """)
    
    return '\n'.join(html_parts)


# ==========================================
# MAIN EXECUTION
# ==========================================

def main():
    print("🔍 Scanning DE analysis results...")
    data = extract_files_by_celltype()
    
    total_contrasts = sum(len(v) for ct_data in data.values() for v in ct_data.values())
    print(f"✅ Found {total_contrasts} contrasts across {len(data)} cell types")
    
    if total_contrasts == 0:
        print(f"⚠️  No DE results found in {DIR_DE}")
        print(f"   Make sure run_DE_analysis.py has completed successfully")
        return
    
    print(f"\n📝 Generating interactive explorer...")
    html_content = generate_html_explorer(data)
    
    # Write to file
    output_path = os.path.join(DIR_DE, OUTPUT_FILE)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Explorer generated: {output_path}")
    print(f"\n🚀 Open in browser: file://{os.path.abspath(output_path)}")


if __name__ == "__main__":
    main()
