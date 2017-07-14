# mpradesigntools
An R package for generating barcoded Massively Parallel Reporter Assay sequences

# Installation

## Dependencies
MPRA Design Tools depends on the Biostrings and BSgenome.Hsapiens.UCSC.hg38 packages from Bioconductor. First install these in R with the following commands:
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
```

The package also makes use of some tidyverse packages which can be installed with the following commands:
```{r}
install.packages(c('dplyr', 'magrittr', 'purrr', 'readr', 'stringr', 'tibble', 'tidyr', 'purrrlyr'))
```

## Package Installation

If you don't have the devtools package installed, install it like so:
```{r}
install.packages("devtools")
```

After that you can install and load MPRA Design Tools with these commands:

```{r}
devtools::install_github('andrewGhazi/mpradesigntools')
library(mpradesigntools)
```

# Use

This is the companion package to the MPRA Design Tools Shiny application available here: https://andrewghazi.shinyapps.io/designmpra/

The Shiny app allows users to interact with MPRA parameters (such as number of barcodes per allele) and see the effect of changing parameters on the assays power. Researchers can use this to decide what parameters best meet their experimental goals.

Currently the main function of MPRA Design Tools package is to design a set of barcoded sequences for MPRA experiments (without overloading our Shiny server!). This is done with the `processVCF` function. It takes roughly 5 seconds + 10ms per barcoded sequence on a relatively modern CPU, so you can estimate the expected job time in seconds as 
```tex
5 + .01 * Number of barcodes per allele * Number of SNPs in VCF * 2 (for ref/alt alleles)
```

## Example
```{r}
processVCF(inputVCF, barcodesPerAllele, upstreamContextRange, downstreamContextRange, fwdPrimer, reversePrimer, filterPatterns = "AATAAA", outPath = <somewhere>/output.tsv)
processVCF(<pathToVCF>, 5, 75, 75, 'ACTGGCCAG', 'CTCGGCGGCC', filterPatterns = "AATAAA", outPath = <somewhere>/output.tsv)
```
