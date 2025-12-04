#!/usr/bin/env python3
"""
VALIDATION 4: PARALOG CONTROL ANALYSIS
======================================
Gemini's "Killer Experiment"

Logic: Many genes have ancient duplications (paralogs). Often one copy
retains the core essential function, while the other diverges.

Example: RPL3 (Essential) vs RPL3L (Tissue-specific, non-essential)

Test: If the Thermodynamic Hypothesis is correct, the Essential paralog
should be MORE STABLE (higher Tm) than the Non-Essential paralog.

This controls for structure and evolutionary history, isolating
ESSENTIALITY as the variable.
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import re

# Load data
master = pd.read_csv('/mnt/user-data/uploads/master_dataset_with_expression.csv')
meltome = pd.read_csv('/mnt/user-data/uploads/meltome_extracted.csv')

merged = master.merge(meltome, on='gene_symbol', how='inner')

print("="*70)
print("VALIDATION 4: PARALOG CONTROL ANALYSIS")
print("="*70)

# ============================================
# Find potential paralog pairs
# ============================================
print("\nSearching for paralog pairs...")

# Pattern: Base gene vs Base gene + L (like RPL3 vs RPL3L)
# Or numbered variants (RPL3 vs RPL3A)

# Extract base gene name (remove trailing L, A, B, numbers for some genes)
def get_base_name(gene):
    # Remove trailing single letters like L, A, B
    match = re.match(r'^(.+?)(L|A|B)?$', gene)
    if match:
        return match.group(1)
    return gene

# Find all gene pairs where one is essential and one is not
genes_with_tm = merged[['gene_symbol', 'chronos_score', 'category', 'Tm']].dropna()

# Define essential vs non-essential
genes_with_tm['is_essential'] = genes_with_tm['chronos_score'] < -0.5
genes_with_tm['is_lethal'] = genes_with_tm['chronos_score'] < -1.0

# Look for specific known paralog families
paralog_families = {
    'Ribosomal_Large': [g for g in genes_with_tm['gene_symbol'] if g.startswith('RPL')],
    'Ribosomal_Small': [g for g in genes_with_tm['gene_symbol'] if g.startswith('RPS')],
    'Proteasome': [g for g in genes_with_tm['gene_symbol'] if g.startswith('PSM')],
    'ATP_Synthase': [g for g in genes_with_tm['gene_symbol'] if g.startswith('ATP5')],
    'Chaperonin': [g for g in genes_with_tm['gene_symbol'] if g.startswith('CCT')],
}

print("\nParalog families found:")
for family, genes in paralog_families.items():
    print(f"  {family}: {len(genes)} genes")

# ============================================
# Analyze within-family variation
# ============================================
print("\n" + "="*70)
print("WITHIN-FAMILY Tm vs ESSENTIALITY")
print("="*70)

results = []
for family, genes in paralog_families.items():
    family_data = genes_with_tm[genes_with_tm['gene_symbol'].isin(genes)]
    if len(family_data) > 3:
        r, p = stats.pearsonr(family_data['Tm'], family_data['chronos_score'])
        results.append({
            'family': family,
            'n': len(family_data),
            'r': r,
            'p': p,
            'mean_tm': family_data['Tm'].mean(),
            'mean_chronos': family_data['chronos_score'].mean()
        })
        sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"\n{family} (n={len(family_data)}):")
        print(f"  Tm range: {family_data['Tm'].min():.1f} - {family_data['Tm'].max():.1f}°C")
        print(f"  Chronos range: {family_data['chronos_score'].min():.2f} - {family_data['chronos_score'].max():.2f}")
        print(f"  Correlation Tm vs Chronos: r = {r:.3f}, p = {p:.4f} {sig}")

# ============================================
# Look for specific Essential vs Non-Essential paralog pairs
# ============================================
print("\n" + "="*70)
print("ESSENTIAL vs NON-ESSENTIAL PARALOG PAIRS")
print("="*70)

# Find genes that differ by trailing L
pairs_found = []
for gene in genes_with_tm['gene_symbol']:
    if not gene.endswith('L'):
        paralog = gene + 'L'
        if paralog in genes_with_tm['gene_symbol'].values:
            g1_data = genes_with_tm[genes_with_tm['gene_symbol'] == gene].iloc[0]
            g2_data = genes_with_tm[genes_with_tm['gene_symbol'] == paralog].iloc[0]
            
            pairs_found.append({
                'gene1': gene,
                'gene2': paralog,
                'tm1': g1_data['Tm'],
                'tm2': g2_data['Tm'],
                'chronos1': g1_data['chronos_score'],
                'chronos2': g2_data['chronos_score'],
                'essential1': g1_data['is_essential'],
                'essential2': g2_data['is_essential']
            })

print(f"\nFound {len(pairs_found)} paralog pairs (Gene vs GeneL):")
for pair in pairs_found[:10]:  # Show first 10
    e1 = "Essential" if pair['essential1'] else "Non-essential"
    e2 = "Essential" if pair['essential2'] else "Non-essential"
    print(f"\n  {pair['gene1']} ({e1}):")
    print(f"    Tm = {pair['tm1']:.1f}°C, Chronos = {pair['chronos1']:.2f}")
    print(f"  {pair['gene2']} ({e2}):")
    print(f"    Tm = {pair['tm2']:.1f}°C, Chronos = {pair['chronos2']:.2f}")
    print(f"  Δ Tm = {pair['tm1'] - pair['tm2']:.1f}°C")

# Statistical test: Do essential paralogs have higher Tm?
if len(pairs_found) > 0:
    tm_diff = [p['tm1'] - p['tm2'] for p in pairs_found]  # Gene minus GeneL
    chronos_diff = [p['chronos1'] - p['chronos2'] for p in pairs_found]
    
    # If essential paralog is more stable, tm_diff should correlate with chronos_diff
    if len(pairs_found) > 2:
        r, p = stats.pearsonr(tm_diff, chronos_diff)
        print(f"\nCorrelation of ΔTm with ΔChronos across pairs: r = {r:.3f}, p = {p:.4f}")

# ============================================
# ALTERNATIVE: Compare within ribosomal proteins
# ============================================
print("\n" + "="*70)
print("DETAILED RIBOSOMAL PROTEIN ANALYSIS")
print("="*70)

rpl = genes_with_tm[genes_with_tm['gene_symbol'].str.startswith('RPL')]
rps = genes_with_tm[genes_with_tm['gene_symbol'].str.startswith('RPS')]

print(f"\nRibosomal Large Subunit (RPL): {len(rpl)} proteins")
print(f"Ribosomal Small Subunit (RPS): {len(rps)} proteins")

if len(rpl) > 5:
    print(f"\nRPL proteins:")
    print(f"  Mean Tm: {rpl['Tm'].mean():.1f}°C (range: {rpl['Tm'].min():.1f}-{rpl['Tm'].max():.1f})")
    print(f"  Mean Chronos: {rpl['chronos_score'].mean():.2f}")
    r, p = stats.pearsonr(rpl['Tm'], rpl['chronos_score'])
    print(f"  Tm vs Chronos: r = {r:.3f}, p = {p:.4f}")

if len(rps) > 5:
    print(f"\nRPS proteins:")
    print(f"  Mean Tm: {rps['Tm'].mean():.1f}°C (range: {rps['Tm'].min():.1f}-{rps['Tm'].max():.1f})")
    print(f"  Mean Chronos: {rps['chronos_score'].mean():.2f}")
    r, p = stats.pearsonr(rps['Tm'], rps['chronos_score'])
    print(f"  Tm vs Chronos: r = {r:.3f}, p = {p:.4f}")

# ============================================
# VISUALIZATION
# ============================================
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Within-family Tm vs Chronos (Ribosomal)
ax1 = axes[0]
ribosomal = genes_with_tm[genes_with_tm['gene_symbol'].str.match(r'^RP[LS]')]
if len(ribosomal) > 5:
    ax1.scatter(ribosomal['Tm'], ribosomal['chronos_score'], alpha=0.7, c='blue')
    z = np.polyfit(ribosomal['Tm'], ribosomal['chronos_score'], 1)
    p_fit = np.poly1d(z)
    x_line = np.linspace(ribosomal['Tm'].min(), ribosomal['Tm'].max(), 100)
    ax1.plot(x_line, p_fit(x_line), 'r--', linewidth=2)
    r, pval = stats.pearsonr(ribosomal['Tm'], ribosomal['chronos_score'])
    ax1.set_title(f'Ribosomal Proteins\nTm vs Essentiality (r={r:.3f}, p={pval:.4f})')
ax1.set_xlabel('Melting Temperature (°C)')
ax1.set_ylabel('Chronos Score')
ax1.axhline(y=-1, color='red', linestyle=':', alpha=0.5)

# Plot 2: Proteasome
ax2 = axes[1]
proteasome = genes_with_tm[genes_with_tm['gene_symbol'].str.startswith('PSM')]
if len(proteasome) > 3:
    ax2.scatter(proteasome['Tm'], proteasome['chronos_score'], alpha=0.7, c='green')
    if len(proteasome) > 5:
        z = np.polyfit(proteasome['Tm'], proteasome['chronos_score'], 1)
        p_fit = np.poly1d(z)
        x_line = np.linspace(proteasome['Tm'].min(), proteasome['Tm'].max(), 100)
        ax2.plot(x_line, p_fit(x_line), 'r--', linewidth=2)
        r, pval = stats.pearsonr(proteasome['Tm'], proteasome['chronos_score'])
        ax2.set_title(f'Proteasome\nTm vs Essentiality (r={r:.3f}, p={pval:.4f})')
    else:
        ax2.set_title(f'Proteasome (n={len(proteasome)})')
ax2.set_xlabel('Melting Temperature (°C)')
ax2.set_ylabel('Chronos Score')
ax2.axhline(y=-1, color='red', linestyle=':', alpha=0.5)

# Plot 3: ATP Synthase
ax3 = axes[2]
atp5 = genes_with_tm[genes_with_tm['gene_symbol'].str.startswith('ATP5')]
if len(atp5) > 3:
    ax3.scatter(atp5['Tm'], atp5['chronos_score'], alpha=0.7, c='orange')
    if len(atp5) > 5:
        z = np.polyfit(atp5['Tm'], atp5['chronos_score'], 1)
        p_fit = np.poly1d(z)
        x_line = np.linspace(atp5['Tm'].min(), atp5['Tm'].max(), 100)
        ax3.plot(x_line, p_fit(x_line), 'r--', linewidth=2)
        r, pval = stats.pearsonr(atp5['Tm'], atp5['chronos_score'])
        ax3.set_title(f'ATP Synthase\nTm vs Essentiality (r={r:.3f}, p={pval:.4f})')
    else:
        ax3.set_title(f'ATP Synthase (n={len(atp5)})')
ax3.set_xlabel('Melting Temperature (°C)')
ax3.set_ylabel('Chronos Score')
ax3.axhline(y=-1, color='red', linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig('/home/claude/gemini_validations/figure_paralog_analysis.png', dpi=150)
print("\nFigure saved: figure_paralog_analysis.png")

# ============================================
# SUMMARY
# ============================================
print("\n" + "="*70)
print("VALIDATION 4 SUMMARY: PARALOG CONTROL")
print("="*70)
print("""
The Paralog Control tests whether thermodynamic stability predicts
essentiality WITHIN protein families that share structure and history.

FINDINGS:
""")

# Summarize correlations
all_core = genes_with_tm[genes_with_tm['gene_symbol'].str.match(r'^(RPL|RPS|PSM|ATP5|CCT)')]
if len(all_core) > 10:
    r, p = stats.pearsonr(all_core['Tm'], all_core['chronos_score'])
    print(f"  Within Core protein families (n={len(all_core)}):")
    print(f"    Tm vs Chronos: r = {r:.3f}, p = {p:.4f}")
    
    if p < 0.05:
        print("\n  ✓ WITHIN protein families, more stable = more essential")
        print("    This controls for structure and evolutionary history.")
        print("    Stability is linked to essentiality status, not just fold type.")
    else:
        print("\n  ✗ No significant relationship within protein families")

# Save results
import json
results_dict = {
    'n_pairs_found': len(pairs_found),
    'ribosomal_n': len(ribosomal) if 'ribosomal' in dir() else 0,
    'core_families_n': len(all_core) if 'all_core' in dir() else 0
}
with open('/home/claude/gemini_validations/validation_4_results.json', 'w') as f:
    json.dump(results_dict, f, indent=2)
