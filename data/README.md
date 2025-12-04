# Data Files

## master_dataset.csv
Complete gene-level dataset with 18,435 human genes.

| Column | Description |
|--------|-------------|
| gene_symbol | HGNC gene symbol |
| chronos_score | Mean CRISPR essentiality score (DepMap 24Q2) |
| category | Internal, External, or Other |
| phylostrata | Evolutionary age (1 = LUCA) |
| network_degree | Protein-protein interaction count |
| median_tpm | Median expression across GTEx tissues |
| n_go_terms | Number of GO annotations |

## core_91_genes.csv
The 91 universally essential genes conserved across all three domains of life.

## meltome_data.csv
Thermal stability (melting temperature) from the Meltome Atlas.

## disorder_data.csv
Intrinsic disorder predictions from IUPred2A.

---

## Data Sources

- **DepMap**: https://depmap.org/portal/
- **Meltome Atlas**: https://meltomeatlas.org/
- **IUPred2A**: https://iupred2a.elte.hu/
- **GTEx**: https://gtexportal.org/
