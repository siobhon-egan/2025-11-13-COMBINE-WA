---
source: Rmd
title: Reproducible analysis
teaching: "10"
exercises: "0"
---

In this exercise we will cover some aspects relevant to generating and most importantly **re-generating** analysis.

We will use the libraries and data we covered in the previous episode [Application for microbiome data](https://siobhonlegan.com/2025-11-13-COMBINE-WA/10-application-microbiome.html).


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

....

## Making reproducible and easy to share reports

Next I have prepared a `.qmd` file that you will download and then "knit" to reproducible a `.html` file. 


:::::::::::::::::::::::::::::::::::::::: callout

### The Situation - Sharing Results

Imagine you need to give your supervisors a summary of your findings early on, instead of spending weeks neatly formatting and exporting plots one by one in R you may like to share this html file early in your analysis to get some general feedback and direction before find tuning your analysis.

::::::::::::::::::::::::::::::::::::::::::::::::::

It is best practice to not re-install the packages each time when making these reports, so ensure you have the required packages installed prior to running.


``` r
# Install your CRAN packages for the qmd
install.packages("gtsummary")
install.packages("reshape2")
install.packages("knitr")
install.packages("dplyr")

## install your Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("microbiome")
BiocManager::install("phyloseq")
```

Download the `.qmd` file


``` r
download.file(
  url = "https://raw.githubusercontent.com/siobhon-egan/2025-11-13-COMBINE-WA/refs/heads/main/episodes/files/microbiome-data-report.qmd",
  destfile = "data/microbiome-report.qmd"
)
```

Have a go at rendering the `.qmd` file.


:::::::::::::::::::::::::::::::::::::::  challenge

### Challenge

Practical opening and "knitting" the `.qmd` file to html output.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Session info

When sharing R Markdown (.Rmd) or Quarto (.qmd) notebooks, it's good practice to include information about your R session at the end of the document. This helps others (and future you!) understand:

- Which version of R was used
- What packages were loaded
- What platform the code was run on



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

- `devtools::session_info()` - more detailed, especially for package dependencies, and identifies where packages were installed from.
- `renv::snapshot()` – for managing project-specific environments.


## Recommendations

Some recommendations for rendering your `.qmd`:

- Hide warnings and messages from your R chunks
- Keeping your code in the output can be good to see what data was used for analysis and if any filtering or transformations were used, however it is not always neccessary and can make it difficult for other "non-coders" to read. In the "YAML" at the top you will see the option `code-fold: true` is used to allow us to toggle the code on and off.
- Include the `sessionInfo()` at the bottom of your file to keep track of package versions used!

---

:::::::::::::::::::::::::::::::::::::::: keypoints

- Saving your data in `.rda` format can make it easier to load and repeat data analysis.
- A notebook in `.rmd` or `.qmd` format and knitted to html can be an efficient way to store and share data workflows and analysis
- Use `sessionInfo()` when using R to generate reports as a handy way to track package versions used to generate analysis and figures.

::::::::::::::::::::::::::::::::::::::::::::::::::
