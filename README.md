# Canids-ExternalLM

This R code is built to replicate the GM analyses of Hebdon et al. (in revision).

The published version of the code will be included in a release and will be noted here. 

## Required Software

The following software is required for reproducing the analysis and figures: 

 - R version 4.4.1 (2024-06-14)
 - RStudio 2024.09.0+375
 - The following updated packages in R available on CRAN: `tidyverse`, `janitor`, `gtools`, `ggplot2`, 
 `patchwork`, `ggrepel`, `geomorph`.
 - The following packages that are available through Github by using the `devtools` package: 
 `SlicerMorphR`, `munchcolors`. Install these by using the following code: 

```
devtools::install_github('SlicerMorph/SlicerMorphR')
devtools::install_github('lindsaywaldrop/munchcolors')
```

## Reproducing the analysis

In RStudio, open the `GM_workbook.RMD` file in the `doc/` folder. Running all lines 
will reproduce the analysis, statistics, and each figure of the manuscript. Much of
the analysis code is located in `src/GM_stats.R`.

## Modification

To modify this code for your own study you will need to replace/change the following:

 1. to update the landmarks include in the analysis you will need to add/replace the landmark files in the "LMs" folder and run a new GPA analysis in SlicerMorph replacing the contents of the Data folder with the output of that analysis.
 2. A classifier sheet updated for all included animals and a column representing the classification schemes of interest. Further modification of the code is required to accomodate new classifciation schemes 
