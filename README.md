# Thermodynamic Constraints on Gene Essentiality

**Thermodynamic Constraints Define a Core Set of 91 Universally Essential Genes Conserved Across All Three Domains of Life**

[![OSF](https://img.shields.io/badge/OSF-Preregistration-green)](https://osf.io/vxyn8)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Repository:** https://github.com/archethura/thermodynamic-gene-essentiality

---

## Overview

This repository contains data, code, and analysis demonstrating that **thermodynamic stability imposes hard constraints on gene essentiality** that precede and constrain adaptive evolution.

### Key Findings

| Finding | Value | p-value |
|---------|-------|---------|
| **Lethality Gap** (Internal vs External) | 1.46 points | < 10⁻¹³¹ |
| **Thermal Stability Gap** | +2.6°C | 0.03 |
| **Structural Order** | 2× less disorder | < 10⁻⁸ |
| **Within-Family Correlation** (RPS) | r = -0.56 | 0.0002 |
| **Partial Correlation** (controlling expression) | r = -0.296 | 0.0002 |
| **LUCA Enrichment** | >500× | < 10⁻⁹⁰ |

### The Core Finding

91 genes are essential across all three domains of life (Bacteria, Archaea, Eukarya). These genes share a distinct **biophysical signature**:

- **Higher thermal stability** (Tm = 54.6°C vs 52.0°C)
- **Lower intrinsic disorder** (13.4% vs 26.2%)
- **Stronger purifying selection** (dN/dS ≈ 0.05 vs 0.25)
- **Longer protein half-life** (143h vs 26h)

**Critical Evidence for Causality:** Within protein families sharing the same fold and evolutionary history, more stable paralogs are more essential (r = -0.56, p = 0.0002). This cannot be explained by structure or phylogeny alone.

---

## Repository Structure

```
thermodynamic-gene-essentiality/
├── README.md                    # This file
├── LICENSE                      # MIT License
├── data/
│   ├── master_dataset.csv       # Gene-level data (18,435 genes)
│   ├── core_91_genes.csv        # The 91 universally essential genes
│   ├── meltome_data.csv         # Thermal stability (Tm) data
│   └── disorder_data.csv        # IUPred2A disorder predictions
├── analysis/
│   ├── 01_lethality_gap.py      # Essentiality distribution analysis
│   ├── 02_partial_correlation.py # Independence from expression
│   ├── 03_paralog_control.py    # Within-family causality test
│   ├── 04_bootstrap.py          # Bootstrap resampling validation
│   └── 05_cross_domain.py       # E. coli, Sulfolobus, Mouse validation
├── figures/
│   └── [manuscript figures]
└── manuscript/
    ├── main_text.docx
    └── supplementary.docx
```

---

## The Statistical Proof

### Independence from Expression (Table 1)

| Analysis Type | Relationship | Control | r | p-value |
|---------------|--------------|---------|---|---------|
| Simple Correlation | Tm vs Essentiality | None | -0.30 | < 0.001 |
| Confound Check 1 | Expression vs Essentiality | None | -0.14 | 0.18 (ns) |
| Confound Check 2 | Tm vs Expression | None | +0.08 | 0.33 (ns) |
| **Partial Correlation** | **Tm vs Essentiality** | **Expression** | **-0.296** | **0.0002** |

The convergence of simple and partial correlations proves thermodynamic stability is **independent of expression level**.

### Paralog Control (The "Killer Experiment")

| Protein Family | n | r | p-value |
|----------------|---|---|---------|
| Ribosomal Small Subunit (RPS) | 40 | -0.56 | 0.0002 |
| Proteasome (PSM) | 48 | -0.36 | 0.011 |
| All Core Families | 142 | -0.25 | 0.003 |

Within families sharing the same fold: **more stable = more essential**.

---

## The Causal Model

```
Thermodynamic Stability → Evolutionary Persistence → Functional Conservation
         ↓                        ↓                          ↓
    High Tm (+2.6°C)      Strong purifying selection    Essential across
    Low disorder (13%)    (dN/dS ≈ 0.05)               3 domains of life
    High pLDDT (+16)      Long half-life (143h)        91 genes at LUCA
```

---

## Data Sources

| Dataset | Source | Link |
|---------|--------|------|
| Gene Essentiality | DepMap CRISPR (24Q2) | [depmap.org](https://depmap.org/portal/) |
| Thermal Stability | Meltome Atlas | [meltomeatlas.org](https://meltomeatlas.org/) |
| Intrinsic Disorder | IUPred2A | [iupred2a.elte.hu](https://iupred2a.elte.hu/) |
| Protein Structure | AlphaFold DB | [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/) |
| E. coli Essentiality | Keio Collection | PMID: 16738554 |

---

## Links

- **Repository:** https://github.com/archethura/thermodynamic-gene-essentiality
- **Pre-registration:** https://osf.io/vxyn8
- **DepMap Data:** https://depmap.org/portal/

---

## Citation

```bibtex
@article{maddox2024thermodynamic,
  title={Thermodynamic Constraints Define a Core Set of 91 Universally 
         Essential Genes Conserved Across All Three Domains of Life},
  author={Maddox, Joseph},
  journal={bioRxiv},
  year={2024}
}
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.
