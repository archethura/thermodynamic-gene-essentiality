#!/usr/bin/env python3
"""
Analysis 1: The Lethality Gap
=============================
Demonstrates the 1.46-point gap between Internal and External genes.

Key Finding:
- Internal genes: mean Chronos = -1.52 (lethal)
- External genes: mean Chronos = -0.06 (neutral)
- Gap: 1.46 points (p < 10^-131, Cohen's d = 2.69)
"""

import pandas as pd
import numpy as np
from scipy import stats

# Load data
df = pd.read_csv('../data/master_dataset.csv')

# Separate categories
internal = df[df['category'] == 'Internal']['chronos_score']
external = df[df['category'] == 'External']['chronos_score']

# Statistics
print("LETHALITY GAP ANALYSIS")
print("=" * 50)
print(f"Internal genes (n={len(internal)}):")
print(f"  Mean Chronos: {internal.mean():.3f}")
print(f"  Median: {internal.median():.3f}")
print(f"  % Lethal (< -1.0): {(internal < -1.0).mean()*100:.1f}%")

print(f"\nExternal genes (n={len(external)}):")
print(f"  Mean Chronos: {external.mean():.3f}")
print(f"  Median: {external.median():.3f}")
print(f"  % Lethal (< -1.0): {(external < -1.0).mean()*100:.1f}%")

# Statistical test
t_stat, p_val = stats.ttest_ind(internal, external)
cohens_d = (internal.mean() - external.mean()) / np.sqrt(
    (internal.var() + external.var()) / 2
)

print(f"\nLETHALITY GAP: {internal.mean() - external.mean():.2f} points")
print(f"t-statistic: {t_stat:.2f}")
print(f"p-value: {p_val:.2e}")
print(f"Cohen's d: {cohens_d:.2f}")
