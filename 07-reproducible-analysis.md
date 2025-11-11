---
source: Rmd
title: Communicating & Reproducing Analysis
teaching: "10"
exercises: "10"
---

In this exercise we will cover some aspects relevant to communicating and most importantly **re**producing analysis.

We will use the libraries and data we covered in the previous episode [Application for microbiome data](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/06-application-microbiome.html).

## Why the emphasis on reproducible analysis

[![Data lifecycle](fig/clipboard-464455724.png)](https://r4ds.had.co.nz/introduction.html)

Research is rarely a straight line. Imagine this common scenario:

-   A research student runs a preliminary analysis and shares results with their supervisor for feedback.
-   The supervisor suggests changes—perhaps a different filtering threshold, an alternative normalization method, or adding new samples.
-   The student now needs to reanalyse the data, often weeks later. Without a clear record of the original steps, this can mean starting from scratch, risking inconsistencies and wasted time.

The same challenge applies to research staff working in teams. When multiple people collaborate on a project, unclear workflows can lead to:

-   Confusion about which version of the analysis is current.
-   Difficulty reproducing figures or tables for publications.
-   Problems when revisiting the project months or years later.

**Reproducibility = Transparency + Efficiency**

Reproducible analysis ensures that:

-   Every step of your workflow is documented.
-   Analyses can be rerun easily after updates or corrections.
-   Collaborators (and future you!) can understand and verify your work.

**Our Approach in This Workshop**

For this exercise, we will use R Markdown and Quarto notebooks—powerful tools for reproducible research. These allow you to:

-   Combine code, results, and narrative in one document.
-   Automatically regenerate outputs when data or code changes.
-   Share complete, executable workflows with colleagues or supervisors.

## Load packages

Now load in the required packages for this analysis


``` r
library(microbiome)
library(knitr)
# load the peerj32 dataset
data(peerj32)
```

## Reproducible data

Let's say you have spent a long time getting your data formatted and read into R.

Now you will want to do some visualizations on this data and it is likely this will need to be repeated (many) times in the future.

You can save your data in and `.rda` or `.rds` format so that it will allow you to quickly read it back into the work space and pick up on more figures. It will also made it easy to read in and generate reports in a reproducible way (which we will do shortly).

Have a go at saving the `peerj32` object that is in your environment to your local directory. If you have followed the previous episodes in this lesson you should already have a **data/** directory.


``` r
save(peerj32, file = "data/peerj32.rda")
```

You can quickly load the `.rda` file back in


``` r
load("data/peerj32.rda")
```

## Making reproducible and easy to share reports

Next you will download a have prepared a `.qmd` file and then work through the code chunks and "knit" the document to a reproducible a `.html` file for communication.

::: callout
### The Situation - Sharing Results

Imagine you need to give your supervisors a summary of your findings early on, instead of spending weeks neatly formatting and exporting plots one by one in R you may like to share this html file early in your analysis to get some general feedback and direction before find tuning your analysis.
:::

It is best practice to not re-install the packages each time when making these reports, so ensure you have the required packages installed prior to running.


``` r
# Install your CRAN packages
install.packages(c(
  "BiocManager",
  "remotes",
  "knitr",
  "dplyr",
  "ggplot2",
  "DT",
  "ggpubr"
))

# Install your Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(
  c("microbiome",
    "phyloseq"
  )
)
```

Download the `.qmd` file using the R code or go [here](https://github.com/siobhon-egan/2025-11-13-COMBINE-WA-afternoon/blob/main/episodes/files/microbiome-data-report.qmd) click the "Download raw file" button.


``` r
# create a directory called files
dir.create("files", showWarnings = FALSE)
# load the qmd into this directory
download.file(
  url = "https://raw.githubusercontent.com/siobhon-egan/2025-11-13-COMBINE-WA/refs/heads/main/episodes/files/microbiome-data-report.qmd",
  destfile = "files/microbiome-report.qmd"
)
```

Have a go at rendering the `.qmd` file.

::: challenge
### Challenge

Practical opening and "knitting" the `.qmd` file to html output.
:::


>You view the rendered html by clicking [here](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/files/microbiome-data-report.html)

::: callout
### When the code looks clean!

In the example on the html document in the heading **Plot taxa prevalence** it is a very simple line of code from the [*Microbiome*](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html) package that generates this figure. The `plot_taxa_prevalence()` does the heavy lifting largely due to the uniform data class structure of the phyloseq object
:::

## Improving your report aesthetics

There are simple strategies you can apply to make your reports more professional, easier to understand, and fully reproducible. These tips apply to both R Markdown (`.Rmd`) and Quarto (`.qmd`) notebooks.

**1. Hide the code**

If you want readers to focus on results rather than code (you may want to also hide warning messages too):

-   *R Markdown*

Add in YAML header:

```         
output: 
  html_document:    
    code_folding: hide
```

-   *Quarto*:

Add in YAML header:

```         
code-fold: true
```

**2. Set global parameters for rendering figures**

Define consistent themes and figure sizes globally: In a setup chunk

```r
knitr::opts_chunk$set(fig.width = 7, fig.height = 5) # set figure sizes
library(ggplot2) # load libraru
theme_set(theme_bw())  # Apply theme to all ggplot
```

**3. Make your code easy to read**

Use RStudio’s built-in code formatting. This automatically indents and spaces your code for better readability. Highlight the code you want to format and press:

-   `CONTROL + SHIFT + A` for windows
-   `COMMAND + SHIFT + A` for mac

**4. Use Meaningful Chunk Names**

Instead of `unnamed-chunk-1`, use descriptive names. This makes debugging and navigation easier. 

Add chunks names in the the curly brackets of your R chunks. For example: `{r data-cleaning}`

## Session info

When sharing R Markdown (.Rmd) or Quarto (.qmd) notebooks, it's good practice to include information about your R session at the end of the document. This helps others (and future you!) understand:

-   Which version of R was used
-   What packages were loaded
-   What platform the code was run on


``` r
sessionInfo()
```

``` output
R version 4.5.2 (2025-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] knitr_1.50        microbiome_1.32.0 ggplot2_4.0.0     phyloseq_1.54.0  

loaded via a namespace (and not attached):
 [1] gtable_0.3.6        xfun_0.54           rhdf5_2.54.0       
 [4] Biobase_2.70.0      lattice_0.22-5      rhdf5filters_1.22.0
 [7] vctrs_0.6.5         tools_4.5.2         generics_0.1.4     
[10] biomformat_1.38.0   stats4_4.5.2        parallel_4.5.2     
[13] tibble_3.3.0        cluster_2.1.8.1     pkgconfig_2.0.3    
[16] Matrix_1.7-4        data.table_1.17.8   RColorBrewer_1.1-3 
[19] S7_0.2.0            S4Vectors_0.48.0    lifecycle_1.0.4    
[22] compiler_4.5.2      farver_2.1.2        stringr_1.6.0      
[25] Biostrings_2.78.0   Seqinfo_1.0.0       codetools_0.2-19   
[28] permute_0.9-8       yaml_2.3.10         pillar_1.11.1      
[31] crayon_1.5.3        tidyr_1.3.1         MASS_7.3-65        
[34] vegan_2.7-2         iterators_1.0.14    foreach_1.5.2      
[37] nlme_3.1-168        tidyselect_1.2.1    digest_0.6.37      
[40] Rtsne_0.17          stringi_1.8.7       purrr_1.2.0        
[43] dplyr_1.1.4         reshape2_1.4.4      splines_4.5.2      
[46] ade4_1.7-23         grid_4.5.2          cli_3.6.5          
[49] magrittr_2.0.4      survival_3.8-3      ape_5.8-1          
[52] withr_3.0.2         scales_1.4.0        XVector_0.50.0     
[55] igraph_2.2.1        multtest_2.66.0     evaluate_1.0.5     
[58] IRanges_2.44.0      mgcv_1.9-1          rlang_1.1.6        
[61] Rcpp_1.1.0          glue_1.8.0          BiocManager_1.30.26
[64] renv_1.1.5          BiocGenerics_0.56.0 jsonlite_2.0.0     
[67] R6_2.6.1            Rhdf5lib_1.32.0     plyr_1.8.9         
```

Some other similar handy functions include:

-   `devtools::session_info()` - more detailed, especially for package dependencies, and identifies where packages were installed from.
-   `renv::snapshot()` – for managing project-specific environments.


------------------------------------------------------------------------

::: keypoints
-   Saving your data in `.rda` format can make it easier to load and repeat data analysis.
-   A notebook in `.rmd` or `.qmd` format and knitted to html can be an efficient way to store and share data workflows and analysis
-   Use `sessionInfo()` when using R to generate reports as a handy way to track package versions used to generate analysis and figures.
- Tips for figures: [*ggpubr*](https://rpkgs.datanovia.com/ggpubr/) package for multipanel figures, consider consistent aesthetics, colouring and themes across figures
:::
