# CRABS: Congruent Rate Analysis in Birth-death Scenarios

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7079514.svg)](https://doi.org/10.5281/zenodo.7079514)

A package for investigating and visualizing the congruence class (Louca & Pennell 2020) of a phylogenetic birth-death model. Previously, this package went by the acronym `ACDC`.

## Installation:

The stable version can be installed from CRAN:

    install.packages("CRABS")

The developmental version can be installed from the github repository. The package `devtools` makes this convenient:

    install.packages("devtools")
    library(devtools)
    install_github("afmagee/CRABS")

## Documentation

There is documentation for each function in R (type `?CRABS` in the R-console), as well as a long-form vignette with an [example workflow](https://afmagee.github.io/CRABS/).

### References
* Louca, S., & Pennell, M. W. (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature, 580(7804), 502-505.
* Höhna, S., Kopperud, B. T., & Magee, A. F. (2022). CRABS: Congruent rate analyses in birth–death scenarios. Methods in Ecology and Evolution, 13, 2709– 2718. https://doi.org/10.1111/2041-210X.13997
* Kopperud, B. T., Magee, A. F., & Höhna, S. (2023). Rapidly Changing Speciation and Extinction Rates Can Be Inferred in Spite of Nonidentifiability. Proceedings of the National Academy of Sciences 120 (7): e2208851120. https://doi.org/10.1073/pnas.2208851120.
* Andreoletti J. & Morlon, H. (2023). Exploring congruent diversification histories with flexibility and parsimony. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14240.
