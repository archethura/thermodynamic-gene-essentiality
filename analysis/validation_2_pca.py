#!/usr/bin/env python3
"""
VALIDATION 2: PCA ANALYSIS
==========================
Addressing the "12 Independent Lines of Evidence" Auto-Correlation

Gemini's concern: The 12 metrics (Tm, disorder, pLDDT, dN/dS, age, etc.)
are likely measuring the same underlying latent variable.

Method: PCA to see if they collapse into a single factor.
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load all available data
master = pd.read_csv('/mnt/user-data/uploads/master_dataset_with_expression.csv')
meltome = pd.read_csv('/mnt/user-data/uploads/meltome_extracted.csv')
disorder = pd.read_csv('/mnt/user-data/uploads/disorder_extracted.csv')

print("="*70)
print("PCA ANALYSIS: DO THE 12 METRICS COLLAPSE INTO ONE FACTOR?")
print("="*70)

# We need to map disorder to gene symbols
# The disorder file has protein_name - need to extract gene symbol
# For now, let's use what we have in master + meltome

# Merge available data
merged = master.merge(meltome, on='gene_symbol', how='inner')
merged = merged.dropna()

print(f"\nGenes with complete data: {len(merged)}")
print(f"Available variables: {merged.columns.tolist()}")

# Define the metrics we can use
# From master: chronos_score, phylostrata, network_degree, median_tpm, n_go_terms
# From meltome: Tm

# Create analysis dataframe
analysis_vars = ['chronos_score', 'phylostrata', 'network_degree', 'median_tpm', 'n_go_terms', 'Tm']
analysis_df = merged[analysis_vars].copy()
analysis_df = analysis_df.dropna()

# Log transform expression
analysis_df['log_tpm'] = np.log10(analysis_df['median_tpm'] + 1)
analysis_df = analysis_df.drop('median_tpm', axis=1)

print(f"\nVariables for PCA: {analysis_df.columns.tolist()}")
print(f"Genes: {len(analysis_df)}")

# ============================================
# CORRELATION MATRIX
# ============================================
print("\n" + "="*70)
print("CORRELATION MATRIX")
print("="*70)

corr_matrix = analysis_df.corr()
print(corr_matrix.round(3))

# ============================================
# PCA ANALYSIS
# ============================================
print("\n" + "="*70)
print("PRINCIPAL COMPONENT ANALYSIS")
print("="*70)

# Standardize
scaler = StandardScaler()
scaled_data = scaler.fit_transform(analysis_df)

# PCA
pca = PCA()
pca_result = pca.fit_transform(scaled_data)

# Variance explained
print("\nVariance Explained by Each Component:")
print("-" * 40)
cumulative = 0
for i, var in enumerate(pca.explained_variance_ratio_):
    cumulative += var
    print(f"  PC{i+1}: {var*100:.1f}% (cumulative: {cumulative*100:.1f}%)")

# Component loadings
print("\nComponent Loadings (correlations with original variables):")
print("-" * 60)
loadings = pd.DataFrame(
    pca.components_.T,
    columns=[f'PC{i+1}' for i in range(len(analysis_df.columns))],
    index=analysis_df.columns
)
print(loadings.round(3))

# ============================================
# INTERPRETATION
# ============================================
print("\n" + "="*70)
print("INTERPRETATION")
print("="*70)

pc1_var = pca.explained_variance_ratio_[0]
pc2_var = pca.explained_variance_ratio_[1]

if pc1_var > 0.5:
    print(f"""
PC1 explains {pc1_var*100:.1f}% of variance - ONE DOMINANT FACTOR.

The metrics are NOT independent - they measure the same thing:
  "Thermodynamic Robustness" = High Tm + Old Age + Low Expression Variability

RECOMMENDATION: Instead of "12 independent lines of evidence,"
claim: "A single latent factor (Thermodynamic Robustness) predicts essentiality."
""")
else:
    print(f"""
PC1 explains only {pc1_var*100:.1f}% - multiple independent factors exist.

PC1: {pc1_var*100:.1f}%
PC2: {pc2_var*100:.1f}%

The metrics capture DIFFERENT aspects of gene properties.
The "12 lines of evidence" claim may be defensible.
""")

# ============================================
# VISUALIZATION
# ============================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Scree plot
ax1 = axes[0, 0]
ax1.bar(range(1, len(pca.explained_variance_ratio_)+1), 
        pca.explained_variance_ratio_*100, 
        alpha=0.7, label='Individual')
ax1.plot(range(1, len(pca.explained_variance_ratio_)+1), 
         np.cumsum(pca.explained_variance_ratio_)*100, 
         'ro-', label='Cumulative')
ax1.set_xlabel('Principal Component')
ax1.set_ylabel('Variance Explained (%)')
ax1.set_title('Scree Plot')
ax1.legend()
ax1.axhline(y=80, color='gray', linestyle='--', alpha=0.5)

# Plot 2: Correlation heatmap
ax2 = axes[0, 1]
im = ax2.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1)
ax2.set_xticks(range(len(corr_matrix.columns)))
ax2.set_yticks(range(len(corr_matrix.columns)))
ax2.set_xticklabels(corr_matrix.columns, rotation=45, ha='right')
ax2.set_yticklabels(corr_matrix.columns)
plt.colorbar(im, ax=ax2, label='Correlation')
ax2.set_title('Correlation Matrix')

# Plot 3: PC1 vs PC2 colored by category
ax3 = axes[1, 0]
merged_with_pca = merged.loc[analysis_df.index].copy()
merged_with_pca['PC1'] = pca_result[:, 0]
merged_with_pca['PC2'] = pca_result[:, 1]

colors = {'Internal': 'blue', 'External': 'orange', 'Other': 'gray'}
for cat in ['Internal', 'External', 'Other']:
    subset = merged_with_pca[merged_with_pca['category'] == cat]
    ax3.scatter(subset['PC1'], subset['PC2'], 
                c=colors[cat], label=cat, alpha=0.5, s=10)
ax3.set_xlabel(f'PC1 ({pc1_var*100:.1f}%)')
ax3.set_ylabel(f'PC2 ({pc2_var*100:.1f}%)')
ax3.set_title('Genes in PC Space')
ax3.legend()

# Plot 4: Loadings biplot
ax4 = axes[1, 1]
for i, var in enumerate(analysis_df.columns):
    ax4.arrow(0, 0, loadings.iloc[i, 0]*3, loadings.iloc[i, 1]*3,
              head_width=0.1, head_length=0.05, fc='red', ec='red')
    ax4.text(loadings.iloc[i, 0]*3.2, loadings.iloc[i, 1]*3.2, var, fontsize=9)
ax4.set_xlim(-1.5, 1.5)
ax4.set_ylim(-1.5, 1.5)
ax4.axhline(y=0, color='gray', linestyle='-', linewidth=0.5)
ax4.axvline(x=0, color='gray', linestyle='-', linewidth=0.5)
ax4.set_xlabel('PC1 Loading')
ax4.set_ylabel('PC2 Loading')
ax4.set_title('Variable Loadings')
ax4.set_aspect('equal')

plt.tight_layout()
plt.savefig('/home/claude/gemini_validations/figure_pca_analysis.png', dpi=150)
print("\nFigure saved: figure_pca_analysis.png")

# ============================================
# PC1 as predictor of essentiality
# ============================================
print("\n" + "="*70)
print("PC1 AS PREDICTOR OF ESSENTIALITY")
print("="*70)

# Correlation between PC1 and Chronos
r, p = stats.pearsonr(pca_result[:, 0], analysis_df['chronos_score'])
print(f"\nPC1 vs Chronos: r = {r:.3f}, p = {p:.2e}")

# By category
for cat in ['Internal', 'External', 'Other']:
    mask = merged.loc[analysis_df.index, 'category'] == cat
    if mask.sum() > 5:
        r, p = stats.pearsonr(pca_result[mask, 0], analysis_df.loc[mask, 'chronos_score'])
        print(f"  {cat}: r = {r:.3f}, p = {p:.2e}, n = {mask.sum()}")

# Save results
import json
results = {
    'pc1_variance': float(pc1_var),
    'pc2_variance': float(pc2_var),
    'cumulative_2pc': float(pc1_var + pc2_var),
    'loadings_pc1': dict(zip(analysis_df.columns, pca.components_[0].tolist())),
    'n_genes': len(analysis_df)
}
with open('/home/claude/gemini_validations/validation_2_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to validation_2_results.json")
