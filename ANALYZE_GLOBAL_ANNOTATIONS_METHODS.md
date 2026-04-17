# Methodology: Global Annotation Analysis

This document describes the mathematical formulas and analytical methods implemented in `analyze_global_annotations.py`.

## 1. Marker Gene Identification
Marker genes for each cell type are identified using a **Student's t-test** (via `scanpy.tl.rank_genes_groups`). 
- **Method**: t-test (overestimated variance).
- **Ranking**: Genes are ranked by their test statistic (score), representing the magnitude of differential expression between the target cell type and all other cells.
- **Top Genes**: The top 10 genes with the highest scores are used for overlap analysis.

## 2. Marker Overlap Ratio
To identify potentially redundant cell type annotations, we calculate the intersection of the top 10 markers between every pair of cell types.
- **Formula**:
  $$Overlap Ratio = \frac{|Markers_{A} \cap Markers_{B}|}{10}$$
- **Threshold**: A ratio $\ge 0.5$ (5 out of 10 genes) is considered a significant overlap.

## 3. Cell Composition Normalization (Heatmap)
To visualize changes in rare cell populations alongside abundant ones, we apply **Z-score normalization** across analysis groups for each cell type.
- **Formula**:
  $$Z_{i,j} = \frac{P_{i,j} - \mu_i}{\sigma_i + \epsilon}$$
  Where:
  - $P_{i,j}$ is the percentage of cell type $i$ in group $j$.
  - $\mu_i$ is the mean percentage of cell type $i$ across all groups.
  - $\sigma_i$ is the standard deviation of cell type $i$ across all groups.
  - $\epsilon$ is a constant ($10^{-6}$) to prevent division by zero.

## 4. Relative Mock-Normalized EBV Response
Instead of a simple subtraction, we use **Relative Fold Change** to measure the impact of EBV infection relative to the Mock baseline for each specific group (Disease/Day).
- **Formula**:
  $$Relative Shift = \frac{P_{EBV} - P_{Mock}}{P_{Mock} + \epsilon}$$
- **Interpretation**: A value of 1.0 represents a 100% increase (doubling) in the cell type proportion compared to its Mock state.

## 5. Temporal Dynamics (Range)
The temporal volatility of a cell type within a condition is measured by the absolute range of its proportions across Day 1, Day 7, and Day 15.
- **Formula**:
  $$Range = \max(P_{Days}) - \min(P_{Days})$$
- **High Dynamics**: Filtered for $Range > 10\%$.

## 6. Statistical Significance
Population shifts between groups are validated using a **two-sample Z-test for proportions** (`statsmodels.stats.proportion.proportions_ztest`).
- **Hypothesis**: $H_0: p_1 = p_2$ (proportions are equal).
- **Outputs**: P-values $< 0.05$ are considered statistically significant.
