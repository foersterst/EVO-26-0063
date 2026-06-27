## Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: phylogenetic comparative evidence from scorpion pincers

### Citation

Please make sure to cite the following paper if you use any data from this repository:

**Foerster, S. Í. A.** (2026). Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: Phylogenetic comparative evidence from scorpion pincers. *Evolution*, qpag097. <https://doi.org/10.1093/evolut/qpag097>

### File Description

- `scripts/data-analysis.R` contains the R code to reproduce all results in the manuscript, including post-processing of RevBayes files.
- `data/scorp-complete.csv` contains linear measurements (species means) for carapace length and chela length, width, and height, all in log<sub>10</sub> scale. Some values in this file were imputed.
- `data/scorp-complete-tree.tre` contains the phylogenetic tree associated with the trait data matrix that includes imputed values (`data/scorp-complete.csv`).
- `data/scorp-complete-noimp.csv` and `data/scorp-complete-noimp.tre` are equivalent to the two files above, but restricted to species with complete trait measurements (i.e., no imputed values).
- `data/data-ssd.csv` contains linear measurements (species means) of carapace length and chela dimensions (length, width, and height) for male and female scorpions. This file is used for the repeatability analyses. Some trait values were imputed separately for each sex.
- The files in `rb/data` are used to run the reparameterized OU models for each chela measurement, with and without imputed values. Chela measurements are expressed as residuals from phylogenetic regressions using carapace length as the predictor, with all variables in log<sub>10</sub> scale.
- `rb/OU-repar.Rev` is the RevBayes script used to run the reparameterized OU models and generate the files in `rb/output`.

### Contact

I welcome inquiries and opportunities for future collaboration. Contact details and information about my research are available at <https://foersterst.github.io>

### Zenodo

```         
[![DOI](https://zenodo.org/badge/1235631639.svg)](https://doi.org/10.5281/zenodo.20973871)
```
