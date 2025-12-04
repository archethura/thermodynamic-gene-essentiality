#!/usr/bin/env python3
"""
VALIDATION 1: RESIDUALS ANALYSIS
================================
The Causal Logic Test

Question: Do these genes persist because they are stable, 
          or are they stable because they must persist?

Method: 
1. Regress Essentiality (Chronos) vs Expression (TPM)
2. Calculate residuals (essentiality not explained by expression)
3. Test if Tm predicts these residuals
4. If yes: Stability contributes to essentiality INDEPENDENT of expression

This decouples the "Dosage Balance" counter-argument.
"""

import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Load data
print("Loading data...")
master = pd.read_csv('/mnt/user-data/uploads/master_dataset_with_expression.csv')
meltome = pd.read_csv('/mnt/user-data/uploads/meltome_extracted.csv')

print(f"Master dataset: {len(master)} genes")
print(f"Meltome data: {len(meltome)} genes")

# Merge datasets
merged = master.merge(meltome, on='gene_symbol', how='inner')
print(f"Merged dataset: {len(merged)} genes with both expression and Tm")

# Filter to genes with valid data
merged = merged.dropna(subset=['chronos_score', 'median_tpm', 'Tm'])
merged = merged[merged['median_tpm'] > 0]  # Need positive expression for log
merged['log_tpm'] = np.log10(merged['median_tpm'] + 1)

print(f"After filtering: {len(merged)} genes")

# Separate Internal vs External
internal = merged[merged['category'] == 'Internal']
external = merged[merged['category'] == 'External']
print(f"Internal genes with Tm: {len(internal)}")
print(f"External genes with Tm: {len(external)}")

# ============================================
# STEP 1: Regress Essentiality vs Expression
# ============================================
print("\n" + "="*60)
print("STEP 1: Essentiality ~ Expression")
print("="*60)

X = sm.add_constant(merged['log_tpm'])
y = merged['chronos_score']

model_expr = sm.OLS(y, X).fit()
print(model_expr.summary().tables[1])

# Calculate residuals
merged['chronos_residual'] = model_expr.resid

print(f"\nExpression explains {model_expr.rsquared*100:.2f}% of essentiality variance")

# ============================================
# STEP 2: Do Tm predict the residuals?
# ============================================
print("\n" + "="*60)
print("STEP 2: Residual Essentiality ~ Thermal Stability (Tm)")
print("="*60)

X_tm = sm.add_constant(merged['Tm'])
y_resid = merged['chronos_residual']

model_tm = sm.OLS(y_resid, X_tm).fit()
print(model_tm.summary().tables[1])

print(f"\nTm explains {model_tm.rsquared*100:.2f}% of RESIDUAL essentiality variance")
print(f"Tm coefficient: {model_tm.params['Tm']:.6f}")
print(f"P-value: {model_tm.pvalues['Tm']:.2e}")

# ============================================
# STEP 3: Full multivariate model
# ============================================
print("\n" + "="*60)
print("STEP 3: Full Model - Essentiality ~ Expression + Tm")
print("="*60)

X_full = sm.add_constant(merged[['log_tpm', 'Tm']])
model_full = sm.OLS(merged['chronos_score'], X_full).fit()
print(model_full.summary().tables[1])

print(f"\nFull model R²: {model_full.rsquared*100:.2f}%")
print(f"Expression-only R²: {model_expr.rsquared*100:.2f}%")
print(f"UNIQUE variance explained by Tm: {(model_full.rsquared - model_expr.rsquared)*100:.2f}%")

# ============================================
# STEP 4: Test within Internal Core only
# ============================================
print("\n" + "="*60)
print("STEP 4: Within Internal Core - Does Tm predict essentiality?")
print("="*60)

if len(internal) > 10:
    X_int = sm.add_constant(internal[['log_tpm', 'Tm']])
    model_int = sm.OLS(internal['chronos_score'], X_int).fit()
    print(model_int.summary().tables[1])
    
    # Correlation within Internal
    r_tm, p_tm = stats.pearsonr(internal['Tm'], internal['chronos_score'])
    print(f"\nWithin Internal Core:")
    print(f"  Tm vs Chronos correlation: r = {r_tm:.3f}, p = {p_tm:.4f}")
else:
    print("Not enough Internal genes with Tm data")

# ============================================
# VISUALIZATION
# ============================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Expression vs Essentiality
ax1 = axes[0, 0]
ax1.scatter(merged['log_tpm'], merged['chronos_score'], alpha=0.3, s=5)
ax1.set_xlabel('log10(TPM + 1)')
ax1.set_ylabel('Chronos Score (Essentiality)')
ax1.set_title(f'Expression vs Essentiality\nR² = {model_expr.rsquared:.3f}')
ax1.axhline(y=-1, color='red', linestyle='--', alpha=0.5)

# Plot 2: Tm vs Residual Essentiality
ax2 = axes[0, 1]
ax2.scatter(merged['Tm'], merged['chronos_residual'], alpha=0.3, s=5)
ax2.set_xlabel('Melting Temperature (°C)')
ax2.set_ylabel('Residual Chronos (after expression)')
ax2.set_title(f'Tm vs Residual Essentiality\nR² = {model_tm.rsquared:.3f}, p = {model_tm.pvalues["Tm"]:.2e}')

# Plot 3: Internal vs External - Tm distribution
ax3 = axes[1, 0]
if len(internal) > 0 and len(external) > 0:
    ax3.hist(internal['Tm'].dropna(), bins=30, alpha=0.7, label=f'Internal (n={len(internal)})', density=True)
    ax3.hist(external['Tm'].dropna(), bins=30, alpha=0.7, label=f'External (n={len(external)})', density=True)
    ax3.set_xlabel('Melting Temperature (°C)')
    ax3.set_ylabel('Density')
    ax3.legend()
    
    # T-test
    t_stat, p_val = stats.ttest_ind(internal['Tm'].dropna(), external['Tm'].dropna())
    ax3.set_title(f'Tm Distribution by Category\nt = {t_stat:.2f}, p = {p_val:.4f}')

# Plot 4: Coefficient comparison
ax4 = axes[1, 1]
coefs = ['Expression\n(log TPM)', 'Thermal Stability\n(Tm)']
values = [model_full.params['log_tpm'], model_full.params['Tm']]
errors = [model_full.bse['log_tpm'], model_full.bse['Tm']]
pvals = [model_full.pvalues['log_tpm'], model_full.pvalues['Tm']]

# Standardize for comparison
std_coefs = [
    model_full.params['log_tpm'] * merged['log_tpm'].std(),
    model_full.params['Tm'] * merged['Tm'].std()
]

bars = ax4.bar(coefs, std_coefs, color=['blue', 'orange'])
ax4.set_ylabel('Standardized Coefficient')
ax4.set_title('Relative Contribution to Essentiality\n(Standardized Coefficients)')
ax4.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

# Add significance stars
for i, p in enumerate(pvals):
    if p < 0.001:
        ax4.text(i, std_coefs[i], '***', ha='center', va='bottom')
    elif p < 0.01:
        ax4.text(i, std_coefs[i], '**', ha='center', va='bottom')
    elif p < 0.05:
        ax4.text(i, std_coefs[i], '*', ha='center', va='bottom')

plt.tight_layout()
plt.savefig('/home/claude/gemini_validations/figure_residuals_analysis.png', dpi=150)
print("\nFigure saved: figure_residuals_analysis.png")

# ============================================
# SUMMARY
# ============================================
print("\n" + "="*60)
print("VALIDATION 1 SUMMARY: RESIDUALS ANALYSIS")
print("="*60)
print(f"""
QUESTION: Is thermal stability an independent predictor of essentiality,
          or just a consequence of high expression?

RESULTS:
  - Expression alone explains {model_expr.rsquared*100:.2f}% of essentiality
  - After controlling for expression, Tm explains {model_tm.rsquared*100:.2f}% of RESIDUAL variance
  - In full model, Tm p-value = {model_full.pvalues['Tm']:.2e}
  - Unique variance explained by Tm: {(model_full.rsquared - model_expr.rsquared)*100:.2f}%

INTERPRETATION:
""")

if model_full.pvalues['Tm'] < 0.05:
    print("  ✓ Thermal stability (Tm) is a SIGNIFICANT independent predictor")
    print("    of gene essentiality, even after controlling for expression level.")
    print("  ✓ This SUPPORTS the Thermodynamic Constraint hypothesis.")
    print("  ✓ Stability is NOT merely a consequence of high expression.")
else:
    print("  ✗ Tm is NOT significant after controlling for expression.")
    print("  ✗ The Dosage Balance hypothesis cannot be ruled out.")

# Save results
results = {
    'expression_r2': model_expr.rsquared,
    'tm_residual_r2': model_tm.rsquared,
    'tm_residual_pvalue': model_tm.pvalues['Tm'],
    'full_model_r2': model_full.rsquared,
    'tm_unique_variance': model_full.rsquared - model_expr.rsquared,
    'tm_pvalue_full': model_full.pvalues['Tm'],
    'tm_coefficient': model_full.params['Tm'],
    'n_genes': len(merged)
}

import json
with open('/home/claude/gemini_validations/validation_1_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print("\nResults saved to validation_1_results.json")
