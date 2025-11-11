---
site: sandpaper::sandpaper_site
---

Material for this has been adapted from the from the Carpentries Incubator lesson [The Bioconductor project](https://carpentries-incubator.github.io/bioc-project/). 

This lesson provides an introduction to the Bioconductor project.

A good understanding of the Bioconductor project is the foundation to efficiently use Bioconductor packages for the analysis and visualization of -omics data using R and RStudio.

We download and install R packages from Bioconductor and other repositories to write workflows and perform analyses.
In order to do so, we first identify packages that are available and relevant to our analysis, and we learn from their documentation the best practices to use them as their authors intended it.
For reproducibility, it is also important to identify and track versions of packages used to perform each analysis.

Sometimes, we encounter bugs in packages that we use.
While it is possible to report bugs to the authors and wait for issue to be fixed,
packages hosted on public repositories offer the chance to inspect the code and contribute or propose fixes ourselves.
In addition to being a great opportunity to develop coding skills, community contributors are very often recognized and credited for their contributions!

In this lesson, you will learn:

- To describe the Bioconductor project beyond software packages.
- To navigate the Bioconductor website to find packages for a particular task.
- To install and update Bioconductor packages.
- To open a package vignette and practice running through the examples that they contain.
- To identify standard classes and methods re-used across Bioconductor packages.
- To modify code and contribute to existing Bioconductor packages.
- Best practices to get help from packages developers and peers.

**Additional components for this workshop**

- Apply bioconductor packages to a biological dataset
  - Inspect how multi-data sets are stored
  - Data visualizations
- Generating & **Re-generating** analysis
  - Saving and loading data quickly
  - Reproducible analysis
  - Reports for communicating results



::::::::::::::::::::::::::::::::::::::::::  prereq

## Prerequisites

This workshop has been written is an **intermediate R** audience and as such we have made assumptions that participants are familiar with R. If you'd like to a fresh on some concepts we recommend the following:

- [R for Reproducible Scientific Analysis](https://swcarpentry.github.io/r-novice-gapminder/)
- [Introduction to data analysis with R and Bioconductor](https://carpentries-incubator.github.io/bioc-intro/index.html).
    

::::::::::::::::::::::::::::::::::::::::::::::::::

You will also need to install the following packages that will be used throughout the lesson - do this while the instructor is talking!


```r
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

# Install bioconductor packages
BiocManager::install(
  c(
    "S4Vectors",
    "Biostrings",
    "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg38.masked",
    "GenomicRanges",
    "rtracklayer",
    "biomaRt",
    "microbiome",
    "phyloseq"
  )
)

# Install packages from github
remotes::install_git("https://codeberg.org/hrbrmstr/hrbrthemes.git")
```

---

This lesson has been made on the [lesson-intro-r-bioconductor](https://carpentries-incubator.github.io/bioc-intro/index.html) and built with [The Carpentries Workbench][workbench]. 

[workbench]: https://carpentries.github.io/sandpaper-docs

