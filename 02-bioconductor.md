---
source: Rmd
title: Bioconductor Introduction & Install
teaching: "15"
exercises: "5"
---



::::::::::::::::::::::::::::::::::::::: objectives

- Describe the Bioconductor project globally.
- Gain a global view of the Bioconductor project in the R ecosystem.
- Identify sources of information to watch for future updates about the Bioconductor project.
- Install BiocManager.
- Install Bioconductor packages.

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::: questions

- What does the Bioconductor project comprise?
- How does the Bioconductor project relate to the CRAN repository?
- How can I learn to use Bioconductor packages effectively?
- How do I join and communicate with the Bioconductor community?
- How do I install Bioconductor packages?
- How do I check if newer versions of my installed packages are available?
- How do I update Bioconductor packages?
- How do I find out the name of packages available from the Bioconductor repositories?

::::::::::::::::::::::::::::::::::::::::::::::::::



## What is Bioconductor?

### A brief history of Bioconductor

The Bioconductor project was started in the Fall of 2001, as an initiative for the collaborative creation of extensible software for computational biology and bioinformatics (Gentleman, Carey, Bates, Bolstad, Dettling, Dudoit, Ellis, Gautier, Ge, Gentry, Hornik, Hothorn, Huber, Iacus, Irizarry, Leisch, Li, Maechler, Rossini, Sawitzki, Smith, Smyth, Tierney, Yang, and Zhang, 2004).
From the very start, the stated mission of the project was to develop tools for the statistical analysis and comprehension of large datasets and technological artifacts in rigorously and robustly designed experiments.
Beyond statistical analyses, the interpretation of statistical results is supported by packages providing biological context, visualization, and reproducibility.

Over the years, software packages contributed to the Bioconductor project have reflected the evolution and emergence of several high-throughput technologies, from microarrays to single-cell genomics, through many variations of sequencing experiments (e.g., RNA-seq, ChIP-seq, DNA-seq), analyses (e.g., sequence variation, copy number variation, single nucleotide polymorphisms), and data modalities (e.g., flow cytometry, proteomics, microscopy and image analysis).

Crucially, the project has not only released software packages implementing novel statistical tests and methodologies, but also produced a diverse range of packages types granting access to databases of molecular annotations and experimental datasets.

The Bioconductor project culminates at an annual conference in North America in the summer, while regional conferences offer great opportunities for networking in Europe, Asia, and North America.
The project is committed to promote a diverse and inclusive community, including a [Code of Conduct][bioc-code-of-conduct] enforced by a Code of Conduct committee.

![](fig/bioc-timeline.svg){alt='Timeline of major Bioconductor milestones alongside technological advancements.'}

**Timeline of major Bioconductor milestones alongside technological advancements.**
Above the timeline, the figure marks the first occurence of major events.
Within the timeline, the name of packages providing core infrastructure indicate the release date.
Below the timeline, major technological advancements contextualise the evolution of the Bioconductor project over time.

### A scientific project

The original publication describes the aims and methods of the project at its inception Gentleman, Carey, Bates et al. (2004).

Huber, Carey, Gentleman, Anders, Carlson, Carvalho, Bravo, Davis, Gatto, Girke, Gottardo, Hahne, Hansen, Irizarry, Lawrence, Love, MacDonald, Obenchain, Oles, Pages, Reyes, Shannon, Smyth, Tenenbaum, Waldron, and Morgan (2015) illustrates the progression of the project, including descriptions of core infrastructure and case studies, from the perspective of both users and developers.

Amezquita, Lun, Becht, Carey, Carpp, Geistlinger, Marini, Rue-Albrecht, Risso, Soneson, Waldron, Pages, Smith, Huber, Morgan, Gottardo, and Hicks (2020) reviewed further developments of the project in the wake of single-cell genomics technologies.

Many more publications and book chapters cite the Bioconductor project, with recent example listed on the [Bioconductor website][bioc-publications].

### A package repository

#### Overview and relationship to CRAN

Undoubtedly, software packages are the best-known aspect of the Bioconductor project.
Since its inception in 2001, the repository has grown over time to host thousands of packages.

The Bioconductor project has extended the preexisting CRAN repository -- much larger and general-purpose in scope -- to comprise R packages primarily catering for bioinformatics and computational biology analyses.

:::::::::::::::::::::::::::::::::::::::::  callout

### Going further

The [Discussion][discuss-cran] article of this lesson includes a section discussing the relationship of Bioconductor and CRAN in further details.

::::::::::::::::::::::::::::::::::::::::::::::::::

#### The Bioconductor release cycle

The Bioconductor project also extended the packaging infrastructure of the CRAN repository to better support the deployment and management of packages at the user level (Gentleman, Carey, Bates et al., 2004).
In particular, the Bioconductor projects features a 6-month release cycle (typically around April and October), which sees a snapshot of the current version of all packages in the Bioconductor repository earmarked for a specific version of R.
R itself is released on an annual basis (typically around April), meaning that for each release of R, two compatible releases of Bioconductor packages are available.

As such, Bioconductor package developers are required to always use the version of R that will be associated with the next release of the Bioconductor project.
This means using the development version of R between October and April, and the release version of R between April and October.

Crucially, the strict Bioconductor release cycle prevents users from installing temporally distant versions of packages that were very likely never tested together.
This practice reflects the development cycle of packages of both CRAN and Bioconductor, where contemporaneous packages are regularly tested by automated systems to ensure that the latest software updates in package dependencies do not break downstream packages, or prompts those package maintainers to update their own software as a consequence.

Prior to each Bioconductor release, packages that do not pass the requires suites of automated tests are deprecated and subsequently removed from the repository.
This ensures that each Bioconductor release provides a suite of packages that are mutually compatible, traceable, and guaranteed to function for the associated version of R.

![](fig/bioc-release-cycle.svg){alt='Timeline of release dates for selected Bioconductor and R versions.'}

**Timeline of release dates for selected Bioconductor and R versions.**
The upper section of the timeline indicates versions and approximate release dates for the R project.
The lower section of the timeline indicates versions and release dates for the Bioconductor project.
Source: [Bioconductor][bioc-release-dates].

#### Package types

Packages are broadly divided in four major categories:

- software
- annotation data
- experiment data
- workflows

[Software packages][glossary-software-package] themselves can be subdivided into packages that provide infrastructure (i.e., classes) to store and access data,
and packages that provide methodological tools to process data stored in those data structures.
This separation of structure and analysis is at the core of the Bioconductor project,
encouraging developers of new methodological software packages to thoughtfully re-use existing data containers where possible,
and reducing the cognitive burden imposed on users who can more easily experiment with alternative workflows without the need to learn and convert between different data structures.

[Annotation data packages][glossary-annotation-package] provide self-contained databases of diverse genomic annotations (e.g., gene identifiers, biological pathways).
Different collections of annotation packages can be found in the Bioconductor project.
They are identifiable by their respective naming pattern, and the information that they contain.
For instance, the so-called `OrgDb` packages (e.g., the *[org.Hs.eg.db](https://bioconductor.org/packages/3.22/org.Hs.eg.db)* package) provide information mapping different types of gene identifiers and pathway databases;
the so-called `EnsDb` (e.g., *[EnsDb.Hsapiens.v86](https://bioconductor.org/packages/3.22/EnsDb.Hsapiens.v86)*) packages encapsulate individual versions of the Ensembl annotations in Bioconductor packages;
and the so-called `TxDb` packages (e.g., *[TxDb.Hsapiens.UCSC.hg38.knownGene](https://bioconductor.org/packages/3.22/TxDb.Hsapiens.UCSC.hg38.knownGene)*) encapsulate individual versions UCSC gene annotation tables.

[Experiment data packages][glossary-experiment-package] provide self-contained datasets that are often used by software package developers to demonstrate the use of their package on well-known standard datasets in their [package vignettes][glossary-vignette].

Finally, [workflow packages][glossary-workflow-package] exclusively provide collections of vignettes that demonstrate the combined usage of several other packages as a coherent workflow, but do not provide any new source code or functionality themselves.

:::::::::::::::::::::::::::::::::::::::  challenge

### Challenge: The Bioconductor website

The Bioconductor website is accessible at [https://bioconductor.org/](https://bioconductor.org/).

Browse the website to find information answering the following questions:

1. How many packages does the current release of the Bioconductor project include?
2. How many packages of each type does this number include?

:::::::::::::::  solution

### Solution

The following solution includes numbers that were valid at the time of writing (Bioconductor release 3.13);
numbers will inevitably be different for future releases of the Bioconductor project.

1. On the page [https://bioconductor.org/](https://bioconductor.org/), in the section "Install", we can read:

> Discover 2042 software packages available in Bioconductor release 3.13.

2. On the page [https://bioconductor.org/](https://bioconductor.org/), in the section "News", click on the link that reads "Bioconductor Bioc `X.Y` Released" (`X.Y` being the version of the current Bioconductor release when you go through this exercise yourself).
  On the linked page, we can read:

> We are pleased to announce Bioconductor 3.13, consisting of 2042 software packages, 406 experiment data packages, 965 annotation packages, and 29 workflows.
> 
> There are 133 new software packages, 22 new data experiment packages, 7 new annotation packages, 1 new workflow, no new books, and many updates and improvements to existing packages; Bioconductor 3.13 is compatible with R 4.1.0, and is supported on Linux, 32- and 64-bit Windows, and macOS 10.14.6 Mojave or higher. This release will include an updated Bioconductor Docker containers.

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

#### Package classification using biocViews

The Bioconductor project uses [biocViews][glossary-biocviews], a set of terms from a controlled vocabulary, to classify Bioconductor packages and facilitate their discovery by thematic search on the [Bioconductor website][biocviews-site].

Each Bioconductor package is tagged with a small set of terms chosen from the available controlled vocabulary, to describe the type and functionality of the package.
Terms are initially selected by the package authors, and subsequently refined during package review or updates to the controlled vocabulary.

:::::::::::::::::::::::::::::::::::::::  challenge

### Challenge

Visit the listing of all packages on the Bioconductor [biocViews][biocviews-site] web page.
Use the "Autocomplete biocViews search" box in the upper left to filter packages by category and explore the graph of software packages by expanding and contracting individual terms.

1. What biocView terms can be used to identify packages that have been tagged for RNA sequencing analysis? ChIP-seq? Epigenetics? Variant annotation? Proteomics? Single-cell genomics?
2. In the `RNASeq` category, two very popular packages are *[DESeq2](https://bioconductor.org/packages/3.22/DESeq2)* and *[edgeR](https://bioconductor.org/packages/3.22/edgeR)*.
  Which one is more popular in terms of download statistics (i.e., lower rank)?

:::::::::::::::  solution

### Solution

1. `RNAseq`, `ChIPSeq`, `Epigenetics`, `VariantAnnotation`, `Proteomics`, `SingleCell`.
2. For Bioconductor release `3.14`, *[DESeq2](https://bioconductor.org/packages/3.22/DESeq2)* and *[edgeR](https://bioconductor.org/packages/3.22/edgeR)* are listed at ranks 23 and 28 respectively.
  In other words, the two packages are among the most frequently downloaded packages in the Bioconductor project, in this instance with a small advantage in favour of *[edgeR](https://bioconductor.org/packages/3.22/edgeR)*.
  
  

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::::::  callout

### Going further

The Bioconductor package *[biocViews](https://bioconductor.org/packages/3.22/biocViews)* is used to support and manage the infrastructure of the controlled vocabulary.
It can also be used to programmatically inspect and subset the list of terms available using their relationship as a graph.

Furthermore, the *[BiocPkgTools](https://bioconductor.org/packages/3.22/BiocPkgTools)* package can be used to browse packages under different [biocViews][glossary-biocviews] (Su, Carey, Shepherd, Ritchie, Morgan, and Davis, 2019).


::::::::::::::::::::::::::::::::::::::::::::::::::

#### Packages interoperability

At the core of the Bioconductor philosophy is the notion of interoperability.
Namely, the capacity of packages to operate on the same data structures.
Importantly, interoperability benefits both users and developers.

Users can more easily write arbitrarily complex workflows that combine multiple packages.
With packages operating on the same data structure, users can maximize their attention to the practical steps of their workflow, and minimize time spent in often complex and error-prone conversions between different data structures specific to each package.
Comparative benchmarks are also easier to implement and can be evaluated more fairly
when competing software packages offering similar functionality operate on input and outputs stored in the same data structures.

Similarly, developers of new packages can focus on the implementation of novel functionality, borrowing existing data structures that offer robust and trusted infrastructure for storage, verification, and indexing of information.

Ultimately, the figure below illustrates how many different Bioconductor packages - as well as base R packages - can be combined to carry out a diverse range of analyses, from importing sequencing data into an R session, to the annotation, integration and visualization of data and results.

![](fig/bioc-sequencing-ecosystem.svg){alt='Sequencing Ecosystem.'}

**Sequencing Ecosystem**
Major data processing steps (blue) and relevant software packages (pink) are listed in the context  of archetypal workflows for various types of genomics analyses.
The sequential relation of workflow steps and software package illustrates the importance of interoperability between software package in order to assemble complete end-to-end workflows.

### Conferences, courses and workshops

The Bioconductor community regularly organizes a number of events throughout the year and across the world. For example:

- The annual BioC summer conference in North America
- Regional conference in winter (e.g. BioC Europe, BioC Asia)
- Summer schools (e.g., CSAMA)
- Online meetings open to all community members (e.g., Bioconductor Developers Forum)

Course materials are regularly uploaded on the [Bioconductor website][bioc-course-materials] following each of those events.
In particular, [online books][bioc-books] are being developed and maintained by community members.

The Bioconductor [YouTube channel][bioc-youtube] is used to publish video recordings of conference presentations including talks and workshops, as well as editions of the regular Bioconductor developers forum (link needed).


### Online communication channels

#### Support site

The Bioconductor [support site][bioc-support-site] provides a platform where users and developers can communicate freely (following the Bioconductor [Code of Conduct][bioc-code-of-conduct]) to discuss issues on a range of subjects, ranging from packages to conceptual questions about best practices.

#### Slack workspace

The Bioconductor [Slack workspace][bioc-slack] is an open space that all community members are welcome to join (for free) and use for rapid interactions.
Currently, the "Pro" pricing plan kindly supported by core funding provides:

- Unlimited message archive
- Unlimited apps
- Group video calls with screen sharing
- Work securely with other organizations using Slack Connect

A wide range of channels have been created to discuss a range of subjects, and community members can freely join the discussion on those channels of create new ones to discuss new subjects.

Important announcements are posted on the `#general` channel.

:::::::::::::::::::::::::::::::::::::::::  callout

### Note

Users are encouraged to use the Bioconductor [support site][bioc-support-site] to raise issues that are relevant to the wider community.
The Slack workspace is often most useful for live discussions, and widely subscribed channels (e.g. `#general`) should be used with moderation.

::::::::::::::::::::::::::::::::::::::::::::::::::

#### Developer Mailing List

The [bioc-devel@r-project.org](mailto:bioc-devel@r-project.org) mailing list is used for communication between package developers, and announcements from the Biocondutor core team.

### A scientific and technical community

- [Scientific Advisory Board (SAB)][bioc-scientific-advisory-board] Meet Annually, External and Internal leader in the field who act as project advisors. No Term limits.
- [Technical Advisory Board (TAB)][bioc-technical-advisory-board]. Meet monthly to consider technical aspects of core infastructure and scientific direction of the project. 15 members, 3 year term. Annual open-to-all elections to rotate members. Current officers are Vince Carey (chair), Levi Waldron (vice Chair) Charlotte Soneson (Secretary).
- [Community Advisory Board (CAB)][bioc-community-advisory-board] Meet monthly to consider community outreach, events, education and training. 15 members, 3 year term. Annual open-to-all elections to rotate members. Current officers are Aedin Culhane (chair), Matt Ritchie (co Chair), Lori Kern (Secretary).
- [Code of Conduct committee][bioc-code-of-conduct]

:::::::::::::::::::::::::::::::::::::::::  callout

### Note

At least 1 member of TAB/CAB sits on both to act at the liason to ensure communication of the board.

::::::::::::::::::::::::::::::::::::::::::::::::::

---

## BiocManager

The *[BiocManager](https://bioconductor.org/packages/3.22/BiocManager)* package is the entry point into the Bioconductor package repository.
Technically, this is the only Bioconductor package distributed on the CRAN repository.

It provides functions to safely install Bioconductor packages and check for available updates.

Once the package is installed, the function `BiocManager::install()` can be used to install packages from the Bioconductor repository.
The function is also capable of installing packages from other repositories (e.g., CRAN), if those packages are not found in the Bioconductor repository first.

![](fig/bioc-install.svg){alt='The package BiocManager is available from the CRAN repository and used to install packages from the Bioconductor repository.'}

**The package BiocManager is available from the CRAN repository and used to install packages from the Bioconductor repository.**
The function `install.packages()` from the base R package `utils` can be used to install the *[BiocManager](https://bioconductor.org/packages/3.22/BiocManager)* package distributed on the CRAN repository.
In turn, the function `BiocManager::install()` can be used to install packages available on the Bioconductor repository.
Notably, the `BiocManager::install()` function will fall back on the CRAN repository if a package cannot be found in the Bioconductor repository.

Install the package using the code below.


``` r
install.packages("BiocManager")
```

:::::::::::::::::::::::::::::::::::::::::  callout

### Going further

A number of packages that are not part of the base R installation also provide functions to install packages from various repositories.
For instance:

- `devtools::install()`
- `remotes::install_bioc()`
- `remotes::install_bitbucket()`
- `remotes::install_cran()`
- `remotes::install_dev()`
- `remotes::install_github()`
- `remotes::install_gitlab()`
- `remotes::install_git()`
- `remotes::install_local()`
- `remotes::install_svn()`
- `remotes::install_url()`
- `renv::install()`

Those functions are beyond the scope of this lesson, and should be used with caution and adequate knowledge of their specific behaviors.
The general recommendation is to use `BiocManager::install()` over any other installation mechanism because it ensures proper versioning of Bioconductor packages.


::::::::::::::::::::::::::::::::::::::::::::::::::

## Bioconductor releases and current version

Once the *[BiocManager](https://bioconductor.org/packages/3.22/BiocManager)* package is installed, the `BiocManager::version()` function displays the version (i.e., release) of the Bioconductor project that is currently active in the R session.


``` r
BiocManager::version()
```

``` output
[1] '3.22'
```

Using the correct version of R and Bioconductor packages is a key aspect of reproducibility.
The *[BiocManager](https://bioconductor.org/packages/3.22/BiocManager)* packages uses the version of R running in the current session to determine the version of Biocondutor packages that can be installed in the current R library.

The Bioconductor project produces two releases each year, one around April and another one around October.
The April release of Bioconductor coincides with the annual release of R.
The October release of Bioconductor continues to use the same version of R for that annual cycle (i.e., until the next release, in April).

![](fig/bioc-release-cycle.svg){alt='Timeline of release dates for selected Bioconductor and R versions.'}

**Timeline of release dates for selected Bioconductor and R versions.**
The upper section of the timeline indicates versions and approximate release dates for the R project.
The lower section of the timeline indicates versions and release dates for the Bioconductor project.
Source: [Bioconductor][bioc-release-dates].

During each 6-month cycle of package development, Bioconductor tests packages for compatibility with the version of R that will be available for the next release cycle.
Then, each time a new Bioconductor release is produced, the version of every package in the Bioconductor repository is incremented, including the package *[BiocVersion](https://bioconductor.org/packages/3.22/BiocVersion)* which determines the version of the Bioconductor project.


``` r
packageVersion("BiocVersion")
```

``` output
[1] '3.22.0'
```

This is the case for every package, even those which have not been updated at all since the previous release.
That new version of each package is earmarked for the corresponding version of R;
in other words, that version of the package can only be installed and accessed in an R session that uses the correct version of R.
This version increment is essential to associate a each version of a Bioconductor package with a unique release of the Bioconductor project.

Following the April release, this means that users must install the new version of R to access the newly released versions of Bioconductor packages.

Instead, in October, users can continue to use the same version of R to access the newly released version of Bioconductor packages.
However, to update an R library from the April release to the October release of Bioconductor, users need to call the function `BiocManager::install()` specifying the correct version of Bioconductor as the `version` option, for instance:


``` r
BiocManager::install(version = "3.14")
```

This needs to be done only once, as the *[BiocVersion](https://bioconductor.org/packages/3.22/BiocVersion)* package will be updated to the corresponding version, indicating the version of Bioconductor in use in this R library.

:::::::::::::::::::::::::::::::::::::::::  callout

### Going further

The [Discussion][discuss-release-cycle] article of this lesson includes a section discussing the release cycle of the Bioconductor project.

::::::::::::::::::::::::::::::::::::::::::::::::::

## Check for updates

The `BiocManager::valid()` function inspects the version of packages currently installed in the user library, and checks whether a new version is available for any of them on the Bioconductor repository.

If everything is up-to-date, the function will simply print `TRUE`.


``` r
BiocManager::valid()
```

``` warning
Warning: 2 packages out-of-date; 0 packages too new
```

``` output

* sessionInfo()

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
[1] bibtex_0.5.1     RefManageR_1.4.0 BiocStyle_2.38.0

loaded via a namespace (and not attached):
 [1] httr_1.4.7          cli_3.6.5           knitr_1.50         
 [4] rlang_1.1.6         xfun_0.54           stringi_1.8.7      
 [7] renv_1.1.5          generics_0.1.4      jsonlite_2.0.0     
[10] glue_1.8.0          backports_1.5.0     plyr_1.8.9         
[13] htmltools_0.5.8.1   rmarkdown_2.29      evaluate_1.0.5     
[16] fastmap_1.2.0       yaml_2.3.10         lifecycle_1.0.4    
[19] BiocManager_1.30.26 stringr_1.6.0       compiler_4.5.2     
[22] Rcpp_1.1.0          timechange_0.3.0    digest_0.6.37      
[25] R6_2.6.1            magrittr_2.0.4      tools_4.5.2        
[28] lubridate_1.9.4     xml2_1.4.1         

Bioconductor version '3.22'

  * 2 packages out-of-date
  * 0 packages too new

create a valid installation with

  BiocManager::install(c(
    "digest", "rmarkdown"
  ), update = TRUE, ask = FALSE, force = TRUE)

more details: BiocManager::valid()$too_new, BiocManager::valid()$out_of_date
```

Conveniently, if any package can be updated, the function generates and displays the command needed to update those packages.
Users simply need to copy-paste and run that command in their R console.

:::::::::::::::::::::::::::::::::::::::::  callout

### Example of out-of-date package library

In the example below, the `BiocManager::valid()` function did not return `TRUE`.
Instead, it includes information about the active user session, and displays the exact call to `BiocManager::install()` that the user should run to replace all the outdated packages detected in the user library with the latest version available in CRAN or Bioconductor.

```
> BiocManager::valid()

* sessionInfo()

R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Big Sur 11.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
[1] BiocManager_1.30.16 compiler_4.1.0      tools_4.1.0         renv_0.14.0        

Bioconductor version '3.13'

  * 18 packages out-of-date
  * 0 packages too new

create a valid installation with

  BiocManager::install(c(
    "cpp11", "data.table", "digest", "hms", "knitr", "lifecycle", "matrixStats", "mime", "pillar", "RCurl",
    "readr", "remotes", "S4Vectors", "shiny", "shinyWidgets", "tidyr", "tinytex", "XML"
  ), update = TRUE, ask = FALSE)

more details: BiocManager::valid()$too_new, BiocManager::valid()$out_of_date

Warning message:
18 packages out-of-date; 0 packages too new 
```

Specifically, in this example, the message tells the user to run the following command to bring their installation up to date:

```
  BiocManager::install(c(
    "cpp11", "data.table", "digest", "hms", "knitr", "lifecycle", "matrixStats", "mime", "pillar", "RCurl",
    "readr", "remotes", "S4Vectors", "shiny", "shinyWidgets", "tidyr", "tinytex", "XML"
  ), update = TRUE, ask = FALSE)
```

::::::::::::::::::::::::::::::::::::::::::::::::::

## Exploring the package repository

The Bioconductor [biocViews][glossary-biocviews], demonstrated in the earlier episode [Introduction to Bioconductor][crossref-intro-biocviews], are a great way to discover new packages by thematically browsing the hierarchical classification of Bioconductor packages.

In addition, the `BiocManager::available()` function returns the complete list of package names that are can be installed from the Bioconductor and CRAN repositories.
For instance the total number of numbers that could be installed using *[BiocManager](https://bioconductor.org/packages/3.22/BiocManager)*


``` r
length(BiocManager::available())
```

``` output
[1] 26647
```

Specifically, the union of current Bioconductor repositories and other repositories on the search path can be displayed as follows.


``` r
BiocManager::repositories()
```

``` output
                                                BioCsoft 
           "https://bioconductor.org/packages/3.22/bioc" 
                                                 BioCann 
"https://bioconductor.org/packages/3.22/data/annotation" 
                                                 BioCexp 
"https://bioconductor.org/packages/3.22/data/experiment" 
                                           BioCworkflows 
      "https://bioconductor.org/packages/3.22/workflows" 
                                               BioCbooks 
          "https://bioconductor.org/packages/3.22/books" 
                                             carpentries 
                    "https://carpentries.r-universe.dev" 
                                     carpentries_archive 
                    "https://carpentries.github.io/drat" 
                                                    CRAN 
                              "https://cran.rstudio.com" 
```

Each repository URL can be accessed in a web browser, displaying the full list of packages available from that repository.
For instance, navigate to [https://bioconductor.org/packages/3.14/bioc](https://bioconductor.org/packages/3.14/bioc).

:::::::::::::::::::::::::::::::::::::::::  callout

### Going further

The function `BiocManager::repositories()` can be combined with the base function `available.packages()` to query packages available specifically from any package repository, e.g. the Bioconductor [software package][glossary-software-package] repository.

```
> db = available.packages(repos = BiocManager::repositories()["BioCsoft"])
> dim(db)
[1] 1948   17
> head(rownames(db))
[1] "a4"          "a4Base"      "a4Classif"   "a4Core"      "a4Preproc"
[6] "a4Reporting"
```

::::::::::::::::::::::::::::::::::::::::::::::::::

Conveniently, `BiocManager::available()` includes a `pattern=` argument, particularly useful to navigate annotation resources (the original use case motivating it).
For instance, a range of [Annotation data packages][glossary-annotation-package] available for the mouse model organism can be listed as follows.


``` r
BiocManager::available(pattern = "*Mmusculus")
```

``` output
 [1] "BSgenome.Mmusculus.UCSC.mm10"        "BSgenome.Mmusculus.UCSC.mm10.masked"
 [3] "BSgenome.Mmusculus.UCSC.mm39"        "BSgenome.Mmusculus.UCSC.mm8"        
 [5] "BSgenome.Mmusculus.UCSC.mm8.masked"  "BSgenome.Mmusculus.UCSC.mm9"        
 [7] "BSgenome.Mmusculus.UCSC.mm9.masked"  "EnsDb.Mmusculus.v75"                
 [9] "EnsDb.Mmusculus.v79"                 "PWMEnrich.Mmusculus.background"     
[11] "TxDb.Mmusculus.UCSC.mm10.ensGene"    "TxDb.Mmusculus.UCSC.mm10.knownGene" 
[13] "TxDb.Mmusculus.UCSC.mm39.knownGene"  "TxDb.Mmusculus.UCSC.mm39.refGene"   
[15] "TxDb.Mmusculus.UCSC.mm9.knownGene"  
```

## Installing packages

The `BiocManager::install()` function is used to install or update packages.

The function takes a character vector of package names, and attempts to install them from the Bioconductor repository.


``` r
BiocManager::install(c("S4Vectors", "BiocGenerics"))
```

However, if any package cannot be found in the Bioconductor repository, the function also searches for those packages in repositories listed in the global option `repos`.

:::::::::::::::::::::::::::::::::::::::::  callout

### Contribute !

Add an example of non-Bioconductor package that can be installed using BioManager.
Preferably, a package that will be used later in this lesson.


::::::::::::::::::::::::::::::::::::::::::::::::::

## Uninstalling packages

Bioconductor packages can be removed from the R library like any other R package, using the base R function `remove.packages()`.
In essence, this function simply removes installed packages and updates index information as necessary.
As a result, it will not be possible to attach the package to a session or browse the documentation of that package anymore.


``` r
remove.packages("S4Vectors")
```

## Getting help

For time we will not go over the full episode on [Getting Help](https://carpentries-incubator.github.io/bioc-project/04-getting-help.html), instead below is a some summary of key points and sources:

**The Bioconductor website**

The main [Bioconductor website][bioc-website] provides a host of resources, all freely available without even the need to install R or any Bioconductor package.

In particular, the [biocViews][biocviews-site] page is a great way to thematically explore the collection of packages and identify Bioconductor packages providing a certain functionality.

Furthermore, the website also collects materials from [courses and conferences][bioc-courses-conferences-materials], including presentations, video recordings, and teaching materials.

**Package landing pages**

Each package accepted in the Bioconductor project is granted a landing page on the main Bioconductor website, e.g. *[S4Vectors](https://bioconductor.org/packages/3.22/S4Vectors)*.

Package landing pages contains useful information that can be consulted without the need to install the package itself.
This is particularly useful while browsing the Bioconductor repository is search of packages suitable for a specific task.

**Package vignettes**

When Bioconductor packages are installed in the user's R library, the vignettes associated with that particular version of the package are also installed on the user's computer.
Those locally installed vignettes are the gold standard reference for the version of the package that is currently installed in the R library and used in the R session.
They can be accessed using the function `browseVignettes()`, for instance:


``` r
browseVignettes("BiocManager")
```

**Package help pages**

Bioconductor requires every user-facing package function to be documented in one of the package help pages (often referred to as "man pages", after the name of the package sub-directory where they are stored).

Help pages can be accessed using the `help()` function or the question mark symbol `?`.


``` r
help(topic = "install", package = "BiocManager")
?BiocManager::install
```

**The Bioconductor support site**

The [Bioconductor support site][bioc-support-site] provides a platform for the community of users and developers to ask questions, and help each other through the doubts and challenges related to Bioconductor packages and analytical workflows.

**Workflow packages**

Bioconductor workflow packages are special in the way that they are only expected to contain vignettes, without any additional code or functionality of their own.
Instead, the vignettes of workflow packages exclusively import functionality from other packages, and demonstrate how to combine functions from those packages into an integrated workflow that users are likely to face in their day-to-day work.

Workflow packages can be browsed in a dedicated section of the [biocViews page][biocviews-workflow].

**Slack workspace**

The [Bioconductor Slack workspace][bioc-slack] was created in 2016.

The workspace can be freely joined using this [Heroku app][bioc-slack-heroku] to generate invitations for individual email addresses.

### How to efficiently ask for help

Most community members provide help voluntarily, on their own spare time and without any form of compensation.
As such, it is important to ask questions as clearly as possible, providing relevant information, to help those volunteers identify the source of the issue as rapidly as possible, thus saving time for both the helper and the original poster.

Depending on the question, some key information includes:

- Operating system
- Version of R
- Version of Bioconductor
- Version of individual packages installed in the R library
- Version of individual packages attached to the current session
- Third-party software, libraries, and compilers installed on the user's system
- Source of packages installed in the R library (e.g., Bioconductor, CRAN, GitHub)
- Code executed in the R session leading up to the issue
- Global options active in the session (accessible using `options()`)

## References

[1] R. A. Amezquita, A. T. L. Lun, E. Becht, et al. "Orchestrating
single-cell analysis with Bioconductor". In: _Nat Methods_ 17.2 (2020),
pp. 137-145. ISSN: 1548-7105 (Electronic) 1548-7091 (Linking). DOI:
10.1038/s41592-019-0654-x.
<https://www.ncbi.nlm.nih.gov/pubmed/31792435>.

[2] R. C. Gentleman, V. J. Carey, D. M. Bates, et al. "Bioconductor:
open software development for computational biology and
bioinformatics". In: _Genome Biol_ 5.10 (2004), p. R80. ISSN: 1474-760X
(Electronic) 1474-7596 (Linking). DOI: 10.1186/gb-2004-5-10-r80.
<https://www.ncbi.nlm.nih.gov/pubmed/15461798>.

[3] W. Huber, V. J. Carey, R. Gentleman, et al. "Orchestrating
high-throughput genomic analysis with Bioconductor". In: _Nat Methods_
12.2 (2015), pp. 115-21. ISSN: 1548-7105 (Electronic) 1548-7091
(Linking). DOI: 10.1038/nmeth.3252.
<https://www.ncbi.nlm.nih.gov/pubmed/25633503>.

[4] S. Su, V. Carey, L. Shepherd, et al. "BiocPkgTools: Toolkit for
mining the Bioconductor package ecosystem [version 1; peer review: 2
approved, 1 approved with reservations] ". In: _F1000Research_ 8.752
(2019). DOI: 10.12688/f1000research.19410.1.

[bioc-code-of-conduct]: https://www.bioconductor.org/about/code-of-conduct/
[bioc-publications]: https://www.bioconductor.org/help/publications/
[discuss-cran]: discuss.html#bioconductor-and-CRAN
[bioc-release-dates]: https://bioconductor.org/about/release-announcements/
[glossary-software-package]: reference.html#software-package
[glossary-annotation-package]: reference.html#annotationdata-package
[glossary-experiment-package]: reference.html#experimentdata-package
[glossary-vignette]: reference.html#vignette
[glossary-workflow-package]: reference.html#workflow-package
[glossary-biocviews]: reference.html#biocviews
[biocviews-site]: https://www.bioconductor.org/packages/release/BiocViews.html
[bioc-course-materials]: https://bioconductor.org/help/course-materials/
[bioc-books]: https://www.bioconductor.org/books/release/
[bioc-youtube]: https://www.youtube.com/user/bioconductor
[bioc-support-site]: https://support.bioconductor.org/
[bioc-slack]: https://bioc-community.herokuapp.com/
[bioc-scientific-advisory-board]: https://bioconductor.org/about/scientific-advisory-board/
[bioc-technical-advisory-board]: https://bioconductor.org/about/technical-advisory-board/
[bioc-community-advisory-board]: https://bioconductor.org/about/community-advisory-board/
[bioc-release-dates]: https://bioconductor.org/about/release-announcements/
[discuss-release-cycle]: discuss.html#the-bioconductor-release-cycle
[glossary-biocviews]: reference.html#biocviews
[crossref-intro-biocviews]: https://carpentries-incubator.github.io/bioc-project/02-introduction-to-bioconductor/index.html#package-classification-using-biocviews
[glossary-software-package]: reference.html#software-package
[glossary-annotation-package]: reference.html#annotationdata-package

:::::::::::::::::::::::::::::::::::::::: keypoints

- R packages are but one aspect of the Bioconductor project.
- The Bioconductor project extends and complements the CRAN repository.
- Different types of packages provide not only software, but also annotations, experimental data, and demonstrate the use of multiple packages in integrated workflows.
- Interoperability beteen Bioconductor packages facilitates the writing of integrated workflows and minimizes the cognitive burden on users.
- Educational materials from courses and conferences are archived and accessible on the Bioconductor website and YouTube channel.
- Different channels of communication enable community members to converse and help each other, both as users and package developers.
- The Bioconductor project is governed by scientific, technical, and advisory boards, as well as a Code of Conduct committee.
- The BiocManager package is available from the CRAN repository.
- `BiocManager::install()` is used to install and update Bioconductor packages (but also from CRAN and GitHub).
- `BiocManager::valid()` is used to check for available package updates.
- `BiocManager::version()` reports the version of Bioconductor currently installed.
- `BiocManager::install()` can also be used to update an entire R library to a specific version of Bioconductor.
- The `browseVignettes()` function is recommended to access the vignette(s) installed with each package.
- The Bioconductor main website contains general information, package documentation, and course materials.

::::::::::::::::::::::::::::::::::::::::::::::::::


