---
source: Rmd
title: Introduction
teaching: 5
exercises: 5
---



::::::::::::::::::::::::::::::::::::::: objectives

- Ensure that participants are using the correct version of R to reproduce exactly the contents of this lesson.
- Download the example files for this lesson.
- Provide overview of lesson and structure

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- Am I using the correct version of R for this lesson?
- Why does my version of R matter?
- How do I obtain the files that are used in this lesson?

::::::::::::::::::::::::::::::::::::::::::::::::::

## Version of R

This lesson was developed and tested with R version 4.5.2 (2025-10-31).

Take a moment to launch RStudio and verify that you are using R version 4.5.x, with `x` being any patch version, e.g. 4.5.2.


``` r
R.version.string
```

``` output
[1] "R version 4.5.2 (2025-10-31)"
```

This is important because Bioconductor uses the version of R running in the current session to determine the version of Bioconductor packages that can be installed in the R library associated with the current R session.
Using a different version of R while following this lesson may lead to unexpected results.

## Download files

Several episodes in this lesson rely on example files that participants need to download.

Run the code below programmatically create a folder called `data` in the current working directory, and download the lesson files in that folder.


``` r
dir.create("data", showWarnings = FALSE)
download.file(
    url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/master/data/TrimmomaticAdapters/TruSeq3-PE-2.fa", 
    destfile = "data/TruSeq3-PE-2.fa"
)
download.file(
    url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/master/data/ActbGtf/actb.gtf", 
    destfile = "data/actb.gtf"
)
download.file(
    url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/master/data/ActbOrf/actb_orfs.fasta", 
    destfile = "data/actb_orfs.fasta"
)
download.file(
  url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/devel/data/SummarizedExperiment/counts.csv",
  destfile = "data/counts.csv"
)
download.file(
  url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/devel/data/SummarizedExperiment/gene_metadata.csv",
  destfile = "data/gene_metadata.csv"
)
download.file(
  url = "https://raw.githubusercontent.com/Bioconductor/bioconductor-teaching/devel/data/SummarizedExperiment/sample_metadata.csv",
  destfile = "data/sample_metadata.csv"
)
```

:::::::::::::::::::::::::::::::::::::::::  callout

## Note

Ideally, participants might want to create a new [RStudio project][external-rstudio-project] and download the lesson files in a sub-directory of that project.

Using an RStudio project sets the working directory to the root directory of that project.
As a consequence, code is executed relative to that root directory, often avoiding the need for using absolute file paths to import/export data from/to files.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Structure of lesson

The content in this lesson has been heavily influenced from the material in the Carpentries Incubator lesson [The Bioconductor project](https://carpentries-incubator.github.io/bioc-project/). 

The first three lessons have been restructured from episodes with additional material added to provide some practice with storing biological data sets, visualization, reproducible workflow and communicating results.

- **Part A - Bioconductor episodes** - a revised and reduced content from the Carpentries Incubator lesson [The Bioconductor project](https://carpentries-incubator.github.io/bioc-project/). 
  - [Introduce, install and getting help with bioconductor](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/02-bioconductor.html)
  - [Working with biological sequences](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/03-biological-sequences.html)
  - [Reproducible annotations](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/04-annotations.html)
  - [Storing and organsing data using the s4 data-class and an example with `summarizedexperiment`](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/05-s4.html)
- **Part B - [application with biological data](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/06-application-microbiome.html)** - work through an example data set and packages from bioconductor for the storage, structure, and visualization of datasets
- **Part C - [reproducibility & communication](https://siobhonlegan.com/2025-11-13-COMBINE-WA-afternoon/07-reproducible-analysis.html)** - a workflow using principals described in previous lessons to work through some practices for reproducible workflows and strategies for communicating and sharing results. 


Material in Part A (the first four episodes) has been trimmed down, and we direct you to the full lesson material available in [The Bioconductor project](https://carpentries-incubator.github.io/bioc-project/). Note that the [Genomic ranges](https://carpentries-incubator.github.io/bioc-project/07-genomic-ranges.html) episode has been removed from content presented here.

Episodes for parts B and C are longer workflows which we will work through together covering some fundamental concepts alongside lots of hands on exercises.

[external-rstudio-project]: https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects


:::::::::::::::::::::::::::::::::::::::: keypoints

- Participants will only be able to install the version of Bioconductor packages described in this lesson and reproduce their exact outputs if they use the correct version of R.
- The files used in this lesson should be downloaded in a local path that is easily accessible from an R session.

::::::::::::::::::::::::::::::::::::::::::::::::::


