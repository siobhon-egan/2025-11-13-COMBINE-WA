---
title: Setup
---

## Software Setup

This lesson assumes you have R and RStudio installed on your computer.

- [Download and install the latest version of R](https://www.r-project.org/).
- [Download and install RStudio](https://www.rstudio.com/products/rstudio/download/#download). RStudio is an application (an integrated development environment or IDE) that facilitates the use of R and offers a number of nice additional features. You will need the free Desktop version for your computer.

If you don't already have R and RStudio installed, follow the instructions for your operating system below below.
You have to install R before you install RStudio. 

This setup guide is taken from, <https://github.com/datacarpentry/R-ecology-lesson/blob/main/learners/setup.md> and modified for use in this lesson.

## Preparations

Data Carpentry's teaching is hands-on, and to follow this lesson learners must have R and RStudio installed on their computers.
They also need to be able to install a number of R packages, create directories, and download files.

To avoid troubleshooting during the lesson, learners should follow the instructions below to download and install everything beforehand.
If the computer is managed by their organization's IT department they might need help from an IT administrator.

### Install R and RStudio

R and RStudio are two separate pieces of software:

- **R** is a programming language and software used to run code written in R.
- **RStudio** is an integrated development environment (IDE) that makes using R easier. In this course we use RStudio to interact with R.

If you don't already have R and RStudio installed, follow the instructions for your operating system below.
You have to install R before you install RStudio.

::::::: spoiler

## For Windows

- Download R from the [CRAN website](https://cran.r-project.org/bin/windows/base/release.htm).
- Run the `.exe` file that was just downloaded
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
- Under _Installers_ select **Windows Vista 10/11 - RSTUDIO-xxxx.yy.z-zzz.exe** (where x = year, y = month, and z represent version numbers)
- Double click the file to install it
- Once it's installed, open RStudio to make sure it works and you don't get any error messages.

:::::::::::::::

::::::: spoiler

## For MacOS

- Download R from the [CRAN website](https://cran.r-project.org/bin/macosx/)
- Select the `.pkg` file for the latest R version
- Double click on the downloaded file to install R
- It is also a good idea to install [XQuartz](https://www.xquartz.org/) (needed by some packages)
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
- Under _Installers_ select **Mac OS 13+ - RSTUDIO-xxxx.yy.z-zzz.dmg** (where x = year, y = month, and z represent version numbers)
- Double click the file to install RStudio
- Once it's installed, open RStudio to make sure it works and that you don't get any error messages

:::::::::::::::

::::::: spoiler

## For Linux

- Click on your distribution in the [Linux folder of the CRAN website](https://cran.r-project.org/bin/linux/). Linux Mint users should follow instructions for Ubuntu.
- Go through the instructions for your distribution to install R.
- Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
- Select the relevant installer for your Linux system (Ubuntu/Debian or Fedora)
- Double click the file to install RStudio
- Once it's installed, open RStudio to make sure it works, and you don't get any error messages.

:::::::::::::::

### Update R and RStudio

If you already have R and RStudio installed, first check if your R version is up to date:

- When you open RStudio your R version will be printed in the console on the bottom left. Alternatively, you can type `sessionInfo()` into the console. If your R version is 4.0.0 or later, you don't need to update R for this lesson. If your version of R is older than that, download and install the latest version of R from the R project website [for Windows](https://cran.r-project.org/bin/windows/base/), [for MacOS](https://cran.r-project.org/bin/macosx/), or [for Linux](https://cran.r-project.org/bin/linux/)
- It is not necessary to remove old versions of R from your system, but if you wish to do so you can check [How do I uninstall R?](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#How-do-I-UNinstall-R_003f)
- After installing a new version of R, you will have to reinstall all your packages with the new version.
  For Windows, there is a package called `installr` that can help you with upgrading your R version and migrate your package library.
- `rig` will allow you to install and switch easily between R versions using Windows, macOS and Linux, and keep your R installation updated.
  See the [documentation](https://github.com/r-lib/rig) for `rig` for more details.
- To update RStudio to the latest version, open RStudio and click on `Help > Check for Updates`.
  If a new version is available follow the instruction on screen. By default, RStudio will also automatically notify you of new versions every once in a while.

::::::::::::::::::::: callout

The changes introduced by new R versions are usually backwards-compatible.
That is, your old code should still work after updating your R version.
However, if breaking changes happen, it is useful to know that you can have multiple versions of R installed in parallel and that you can switch between them in RStudio by going to `Tools > Global Options > General > Basic`.

While this may sound scary, it is **far more common** to run into issues due to using out-of-date versions of R or R packages.
Keeping up with the latest versions of R, RStudio, and any packages you regularly use is a good practice.

:::::::::::::::::::::::::::::

## Install R packages

You will also need to install the following packages that will be used throughout the lesson - do this while the instructor is talking!
lesson.

```r
# specify latest BiocManager version compatible with R 4.5.1
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.22")

# Install rest of your CRAN packages
install.packages(c(
  "remotes",
  "knitr",
  "dplyr",
  "ggplot2",
  "DT",
  "ggpubr"
))

# Check CRAN library installs
library(BiocManager)
library(remotes)
library(knitr)
library(dplyr)
library(ggplot2)
library(DT)
library(ggpubr)


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

# Check bioconductor library installs
library(S4Vectors)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38.masked)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(microbiome)
library(phyloseq)
```

:::::::::::::::::::::::::::::::::::::::::  callout
If you have issues downloading packages try installing them one at a time.
::::::::::::::::::::::::::::::::::::::::::::::::::
