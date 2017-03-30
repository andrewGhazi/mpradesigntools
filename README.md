# mpradesigntools
A tool for generating barcoded Massively Parallel Reporter Assay sequences

# Installation

MPRA Design Tools depends on the Biostrings and BSgenome.Hsapiens.UCSC.hg38 packages from Bioconductor. First install these with the following commands:
```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
```

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

Currently the main function of MPRA Design Tools is to design a set of barcoded sequences for MPRA experiments. This is done with the `processVCF` function. It takes roughly 

## Example
```{r}
processVCF(<pathToVCF>, 5, 75, 'ACTGGCCAG', 'CTCGGCGGCC', outPath = <somewhere>/output.tsv)
```
