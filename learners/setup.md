---
title: Setup
---

## Software Setup

This lesson assumes you have R and RStudio installed on your computer.

- [Download and install the latest version of R](https://www.r-project.org/).
- [Download and install RStudio](https://www.rstudio.com/products/rstudio/download/#download). RStudio is an application (an integrated development environment or IDE) that facilitates the use of R and offers a number of nice additional features. You will need the free Desktop version for your computer.

If you don't already have R and RStudio installed, follow the instructions for your operating system below.
You have to install R before you install RStudio. 

::::::: spoiler

### For Windows

* Download R from the [CRAN website](https://cran.r-project.org/bin/windows/base/release.htm).
* Run the `.exe` file that was just downloaded
* Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
* Under *Installers* select **Windows Vista 10/11 - RSTUDIO-xxxx.yy.z-zzz.exe** (where x = year, y = month, and z represent version numbers)
* Double click the file to install it
* Once it's installed, open RStudio to make sure it works and you don't get any error messages.
  
:::::::::::::::::::::::::

:::::::::::::::: spoiler

### For MacOS

* Download R from the [CRAN website](https://cran.r-project.org/bin/macosx/).
* Select the `.pkg` file for the latest R version
* Double click on the downloaded file to install R
* It is also a good idea to install [XQuartz](https://www.xquartz.org/) (needed by some packages)
* Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
* Under *Installers* select **Mac OS 13+ - RSTUDIO-xxxx.yy.z-zzz.dmg** (where x = year, y = month, and z represent version numbers)
* Double click the file to install RStudio
* Once it's installed, open RStudio to make sure it works and you don't get any error messages.

::::::::::::::::

::::::: spoiler

### For Linux

* Click on your distribution in the [Linux folder of the CRAN website](https://cran.r-project.org/bin/linux/). Linux Mint users should follow instructions for Ubuntu.
* Go through the instructions for your distribution to install R.
* Go to the [RStudio download page](https://www.rstudio.com/products/rstudio/download/#download)
* Select the relevant installer for your Linux system (Ubuntu/Debian or Fedora)
* Double click the file to install RStudio
* Once it's installed, open RStudio to make sure it works and you don't get any error messages.

::::::::::::::::

## Update R and RStudio

If you already have R and RStudio installed, first check if your R version is up to date:

* When you open RStudio your R version will be printed in the console on the bottom left. Alternatively, you can type `sessionInfo()` into the console. If your R version is 4.0.0 or later, you don't need to update R for this lesson. If your version of R is older than that, download and install the latest version of R from the R project website [for Windows](https://cran.r-project.org/bin/windows/base/), [for MacOS](https://cran.r-project.org/bin/macosx/), or [for Linux](https://cran.r-project.org/bin/linux/)
* It is not necessary to remove old versions of R from your system, but if you wish to do so you can check [How do I uninstall R?](https://cran.r-project.org/bin/windows/base/rw-FAQ.html#How-do-I-UNinstall-R_003f) 
* After installing a new version of R, you will have to reinstall all your packages with the new version. For Windows, there is a package called `installr` that can help you with upgrading your R version and migrate your package library. A similar package called `pacman` can help with updating R packages across
To update RStudio to the latest version, open RStudio and click on 
`Help > Check for Updates`. If a new version is available follow the 
instruction on screen. By default, RStudio will also automatically notify you 
of new versions every once in a while.

### Updating R packages

Generally, it is recommended to keep your R version and all packages up to date, because new versions bring improvements and important bugfixes. To update the packages that you have installed, click `Update` in the `Packages` tab in the bottom right panel of RStudio, or go to `Tools > Check for Package Updates...` 

You should update **all of the packages** required for the lesson, even if you installed them relatively recently.

Sometimes, package updates introduce changes that break your old code, which can be very frustrating. To avoid this problem, you can use a package called `renv`. It locks the package versions you have used for a given project and makes it straightforward to reinstall those exact package version in a new environment, for example after updating your R version or on another computer. However, the details are outside of the scope of this lesson.

