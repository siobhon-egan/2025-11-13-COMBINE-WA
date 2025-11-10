---
source: Rmd
title: Microbiome application
teaching: "10"
exercises: "0"
---

In this example we will explore some bioconductor packages for the analysis of microbiome data.

During this we will be making use of two main packages from the bioconductor library:
- [*Microbiome*](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html)
- [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
-

We will also introduce a new data class - `phyloseq-class` from the package [*Phyloseq*](https://joey711.github.io/phyloseq/) that uses the S4 class system.

> The [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) package is designed to store and manage phylogenetic sequencing data, including abundance data, phylogenetic trees, taxonomic assignments, and sample metadata, within a single experiment-level S4 object. This S4 object-oriented approach provides a structured and organized way to handle complex microbiome data, facilitating reproducible analyses and data sharing.

At its core is the `phyloseq-class` object, which integrates multiple data types: OTU/ASV tables, taxonomic annotations, sample metadata, and phylogenetic trees.
While originally developed for 16S rRNA microbiome studies, phyloseq-class objects are also widely used in:
 
 - eDNA metabarcoding studies (e.g., biodiversity monitoring in aquatic or terrestrial environments)
- Fungal ITS sequencing projects
- Metagenomic and metatranscriptomic datasets with annotated taxonomic profiles
- Environmental genomics and ecological community profiling

## Install packages

Install our [*Microbiome*](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html) and [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) packages using `BiocManager::install()`


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

We will explore some example datasets that come with the packages [Microbiome*](https://microbiome.github.io/tutorials/Data.html) and [*Phyloseq*](https://joey711.github.io/phyloseq/Example-Data). We will take some time to explore what information they contain, how we can view the datasets and how they differ.

**Explore dataset 1 - `dietswap`**

The first dataset we will load is from [O'Keefe S et al. Nature Communications, 2015](https://dx.doi.org/10.1038/ncomms7342) made available by [*Microbiome*](https://microbiome.github.io/tutorials/Data.html) package. This study 
contains microbiome data from a study with African and African American groups undergoing a two-week diet swap. This data set is based on the Human Intestinal Tract (HIT)Chip phylogenetic 16S microarray ([Rajilic-Stojanovic et al. 2009](https://doi.org/10.1111/j.1462-2920.2009.01900.x)). This profiling technology differs from the more widely used 16S rRNA amplicon sequencing. Column metadata includes the subject identifier, sex, nationality, group information, sample identifier, time point information, time point information within group and BMI group. Row metadata contains taxonomic information on the Phylum, Family and Genus level.


``` r
data(dietswap)
```

**Explore dataset 2 - `peerj32`**

The second dataset we will load is the `peerj32` data from [Lahti et al. PeerJ 1:e32, 2013](https://doi.org/10.7717/peerj.32) made available by [*Microbiome*](https://microbiome.github.io/tutorials/Data.html) package. This characterizes associations between human intestinal microbiota and blood serum lipids. Note that this data set contains an additional data matrix of lipid species. 


``` r
data(peerj32)
```

**Explore dataset 3 - `GlobalPatterns`**

The third dataset we will load is the `GlobalPatterns` data from [Caporaso et al. PNAS, 108, 4516-4522, 2011](https://doi.org/10.1073/pnas.1000080107) made available by [*Phyloseq*](https://joey711.github.io/phyloseq/Example-Data). This work compared the microbial communities from 25 environmental samples and three known "mock communities" -- a total of 9 sample types -- at a depth averaging 3.1 million reads per sample. Authors were able to reproduce diversity patterns seen in many other published studies, while also investigating technical issues/bias by applying the same techniques to simulated microbial communities of known composition.



``` r
data(GlobalPatterns)
```

:::::::::::::::::::::::::::::::::::::::  challenge

### Exploring data

Explore the datasets `dietswap`, `peerj32` and `GlobalPatterns` that you just loaded. You can do this by navigating to your "Environment" pane, or in the console try the following `dietswap@` and see what options are available to you. Can you view the metadata for this data set? Why does `peerj32@` not work in the same way? 

:::::::::::::::  solution

### Solution

To see the full data within the `peerj32` object try typing `peerj32$` into the console. 

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::  challenge

### Missing information

Can you identify what information is missing in each of the three `phyloseq-class` objects for the datasets above?

:::::::::::::::  solution

### Solution

Check the contents of each of the three datasets separately and you can see that both `dietswap` and `peerj32$phyloseq` objects only contain the minimum data required - `otu_table()`, `sample_data()` and `tax_table()`.
The `GlobalPatterns` data also contained additional taxonomic ranks and `phy_tree()` which is a phylogenetic tree that we can use when performing analysis. 

```r
> dietswap
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
```

```r
> peerj32$phyloseq
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 130 taxa and 44 samples ]
sample_data() Sample Data:       [ 44 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
```

```r
> GlobalPatterns
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


## Visualizing complex data

...


:::::::::::::::::::::::::::::::::::::::: keypoints

- Bespoke often domain specific data classes are often built to serve a specific set of data, however by utilising fundamental R principals such as S4 class structures transferring formats between packages for analysis is much easier.
- Look for widely accepted data classes in your field before formatting your data for a specific package early on.
- Use `sessionInfo()` when using R to generate reports as a handy way to track package versions used to generate analysis and figures.

::::::::::::::::::::::::::::::::::::::::::::::::::
