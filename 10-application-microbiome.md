---
source: Rmd
title: Microbiome application
teaching: "10"
exercises: "0"
---

In this example we will explore some bioconductor packages for the analysis of microbiome data.

During this we will also introduce a new data class - `phyloseq-class` from the package [*Phyloseq*](https://joey711.github.io/phyloseq/) that uses the S4 class system.

> The [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) package is designed to store and manage phylogenetic sequencing data, including abundance data, phylogenetic trees, taxonomic assignments, and sample metadata, within a single experiment-level S4 object. This S4 object-oriented approach provides a structured and organized way to handle complex microbiome data, facilitating reproducible analyses and data sharing.

## Install packages

Install our [*microbiome*](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html) package using `BiocManager::install()`


``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("phyloseq")
BiocManager::install("microbiome")
```


Now load in the packages we need for this analysis

``` r
library(microbiome)
library(phyloseq)
library(knitr)
```

## Example datasets

Next we will explore two example datasets that come with the microbiome package. We will take some time to explore what information they contain, how we can view the datasets and how they differ.

**Explore dataset 1 - 'dietswap"**



``` r
data(dietswap)
```

**Explore dataset 2 - 'peerj32"**

Next we will load another data set that contains microbiome data. The `peerj32` data set is from [Lahti et al. PeerJ 1:e32, 2013](https://doi.org/10.7717/peerj.32) characterizes associations between human intestinal microbiota and blood serum lipids. Note that this data set contains an additional data matrix of lipid species. 


``` r
data(peerj32)
```


:::::::::::::::::::::::::::::::::::::::  challenge

### Exploring data

Explore the two datasets `dietswap` and `perj32` that you just loaded. You can do this by navigating to your "Environment" pane, or in the console try the following `pseq@` and see what options are available to you. Can you view the metadata for this data set? Why does `peerj32@` not work in the same way?

:::::::::::::::  solution

### Solution

To see the full data within the `peerj32` object try typing `peerj32$` into the console. 

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


## Visualizing complex data

...


:::::::::::::::::::::::::::::::::::::::: keypoints

- Bespoke often domain specific data classes are often built to serve a specific set of data, however by utilising fundamental R principals such as S4 class structures transferring formats between packages for analysis is much easier.
- Look for widely accepted data classes in your field before formatting your data for a specific package early on.
- Use `sessionInfo()` when using R to generate reports as a handy way to track package versions used to generate analysis and figures.

::::::::::::::::::::::::::::::::::::::::::::::::::
