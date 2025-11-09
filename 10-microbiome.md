---
source: Rmd
title: Microbiome application
teaching: "10"
exercises: "0"
---

In this example we will explore some bioconductor packages for the analysis of microbiome data.

During this we will also introduce a new data class - 'phyloseq' that used the S4 class system.

> The [phyloseq](https://joey711.github.io/phyloseq/) package is designed to store and manage phylogenetic sequencing data, including abundance data, phylogenetic trees, taxonomic assignments, and sample metadata, within a single experiment-level S4 object. This S4 object-oriented approach provides a structured and organized way to handle complex microbiome data, facilitating reproducible analyses and data sharing.

## Install packages

Install our `microbiome` package using `BiocManager`


``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')


BiocManager::install("microbiome")
```


Now load in the packages we need for this analysis

``` r
library(microbiome)
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

### Challenge - exploring data

Explore the two datasets `dietswap` and `perj32` that you just loaded. You can do this by navigating to your "Environment" pane, or in the console try the following `pseq@` and see what options are available to you. Can you view the metadata for this data set? Why does `peerj32@` not work in the same way?

:::::::::::::::  solution

### Solution

To see the full data within the `peerj32` object try typing `peerj32$` into the console. 

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


## Visualizing complex data

...


## Saving your data

Let's say you have spent a long time getting your data formatted and read into R.

Now you will want to do some visualizations on this data and it is likely this will need to be repeated (many) times in the future.

You can save your data in and `.rda` or `.rds` format so that it will allow you to quickly read it back into the work space and pick up on more figures. It will also made it easy to read in and generate reports in a reproducible way (which we will do shortly).

Have a go at saving the `peerj32` object that is in your environment to your local directory.


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

Have a go at rendering the `.qmd` file. 

:::::::::::::::  solution

### Solution

Take note of the `sessionInfo()` at the bottom - this is a handy way to track versions of packages used to generate analysis.


:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


Some recommendations for rendering your `.qmd`:

- Hide warnings and messages from your R chunks
- Keeping your code in the output can be good to see what data was used for analysis and if any filtering or transformations were used, however it is not always neccessary and can make it difficult for other "non-coders" to read. In the "YAML" at the top you will see the option `code-fold: true` is used to allow us to toggle the code on and off.
- Include the `sessionInfo()` at the bottom of your file to keep track of package versions used!

---

:::::::::::::::::::::::::::::::::::::::: keypoints

- Bespoke often domain specific data classes are often built to serve a specific set of data, however by utilising fundamental R principals such as S4 class structures transferring formats between packages for analysis is much easier.
- Look for widely accepted data classes in your field before formatting your data for a specific package early on.
- Use `sessionInfo()` when using R to generate reports as a handy way to track package versions used to generate analysis and figures.

::::::::::::::::::::::::::::::::::::::::::::::::::
