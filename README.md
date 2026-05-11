# Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: phylogenetic comparative evidence from scorpion pincers
This repository accompanies the article Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: phylogenetic comparative evidence from scorpion pincers.

### File Description
- `scripts/data-analysis.R` contains the R code to reproduce all results in the manuscript, including post-processing of RevBayes files.
- `data/scorp-complete.csv` contains linear measurements (species means) for carapace length and chela length, width, and height, all in log<sub>10</sub> scale. Some values in this file were imputed.
- `data/scorp-complete-tree.tre` contains the phylogenetic tree associated with the trait data matrix that includes imputed values (`data/scorp-complete.csv`).
- `data/scorp-complete-noimp.csv` and `data/scorp-complete-noimp.tre` are equivalent to the two files above, but restricted to species with complete trait measurements (i.e., no imputed values).
- `data/data-ssd.csv` contains linear measurements (species means) of carapace length and chela dimensions (length, width, and height) for male and female scorpions. This file is used for the repeatability analyses. Some trait values were imputed separately for each sex.
- The files in `rb/data` are used to run the reparameterized OU models for each chela measurement, with and without imputed values. Chela measurements are expressed as residuals from phylogenetic regressions using carapace length as the predictor, with all variables in log<sub>10</sub> scale.
- `rb/OU-repar.Rev` is the RevBayes script used to run the reparameterized OU models and generate the files in `rb/output`.

### Data usage notice

The data contained in this repository are associated with the publication:

"Biomechanical trade-offs and the evolutionary dynamics of weapon morphology: phylogenetic comparative evidence from scorpion pincers (Evolution, manuscript EVO-26-0063)"

The code in this repository is released under the MIT License.
However, the data files are **not** covered by the MIT License.

The data are provided for transparency and reproducibility purposes only.
Users wishing to reuse, redistribute, or publish analyses based on these data
must first obtain permission from the author.

Please contact:
stenioit@gmail.com

