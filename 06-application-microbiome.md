---
source: Rmd
title: Microbiome application
teaching: "15"
exercises: "0"
---

In this example we will explore some bioconductor packages for the analysis of microbiome data. We will focus on managing and storing data in a structures and reproducible way and then fine tune some visualizations. We will follow on with a similar set on analysis in the next session for [communicating & reproducing analysis]()

During this we will be making use of two main packages from the bioconductor library:

- [*Microbiome*](https://www.bioconductor.org/packages/release/bioc/html/microbiome.html)
- [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)

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

# Working with data

Now we will have a go at working through some data **visualizations** and improving aesthetics ready for publications or sharing with colleagues.

For this exercise we will work through using the `GlobalPatterns` data. 

>The `GlobalPatterns` data is from [Caporaso et al. PNAS, 108, 4516-4522, 2011](https://doi.org/10.1073/pnas.1000080107) made available by [*Phyloseq*](https://joey711.github.io/phyloseq/Example-Data). This work compared the microbial communities from 25 environmental samples and three known "mock communities"

First a re-fresh of what the sample data contains for this dataset


``` r
#show just the first five rows
head(GlobalPatterns@sam_data)
```

``` output
Sample Data:        [6 samples by 7 sample variables]:
        X.SampleID  Primer Final_Barcode Barcode_truncated_plus_T
CL3            CL3 ILBC_01        AACGCA                   TGCGTT
CC1            CC1 ILBC_02        AACTCG                   CGAGTT
SV1            SV1 ILBC_03        AACTGT                   ACAGTT
M31Fcsw    M31Fcsw ILBC_04        AAGAGA                   TCTCTT
M11Fcsw    M11Fcsw ILBC_05        AAGCTG                   CAGCTT
M31Plmr    M31Plmr ILBC_07        AATCGT                   ACGATT
        Barcode_full_length SampleType
CL3             CTAGCGTGCGT       Soil
CC1             CATCGACGAGT       Soil
SV1             GTACGCACAGT       Soil
M31Fcsw         TCGACATCTCT      Feces
M11Fcsw         CGACTGCAGCT      Feces
M31Plmr         CGAGTCACGAT       Skin
                                       Description
CL3       Calhoun South Carolina Pine soil, pH 4.9
CC1       Cedar Creek Minnesota, grassland, pH 6.1
SV1     Sevilleta new Mexico, desert scrub, pH 8.3
M31Fcsw    M3, Day 1, fecal swab, whole body study
M11Fcsw   M1, Day 1, fecal swab, whole body study 
M31Plmr    M3, Day 1, right palm, whole body study
```

Next we often need to do some data pruning or subsetting. In this case we are making sure we remove taxa that doesn't have any sequences present.


``` r
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
```

``` warning
Warning in prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns): 'prune_species' is deprecated.
Use 'prune_taxa' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` warning
Warning in speciesSums(GlobalPatterns): 'speciesSums' is deprecated.
Use 'taxa_sums' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` r
# compare the two datasets
GlobalPatterns
```

``` output
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

``` r
GP
```

``` output
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 18988 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
tax_table()   Taxonomy Table:    [ 18988 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 18988 tips and 18987 internal nodes ]
```

## Alpha diversity plots

For microbiome (or community analysis) a common strategy is to visualize the alpha diversity of samples. Essentially this is how many unique taxa were identified in each sample.

At the simplest level we can plot all the alpha diversity measures from the function `plot_richness()` from the [*Phyloseq*](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html) package. 



``` r
plot_richness(GP)
```

``` warning
Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
ℹ Please use tidy evaluation idioms with `aes()`.
ℹ See also `vignette("ggplot2-in-packages")` for more information.
ℹ The deprecated feature was likely used in the phyloseq package.
  Please report the issue at <https://github.com/joey711/phyloseq/issues>.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
generated.
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

Note that in this case, the Fisher calculation results in a warning (but still plots).

There is a lot going on here so let's try just look at two different alpha diversity measures.

We can avoid this by specifying a measures argument to plot_richness, which will include just the alpha-diversity measures that we want.


``` r
plot_richness(GP, measures=c("Chao1", "Shannon"))
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

:::::::::::::::::::::::::::::::::::::::::  callout

### Getting help on functions

You can find more information on the funciton by typing `?plot_richness` into the console. Have a look and see what measures are available using this function.

::::::::::::::::::::::::::::::::::::::::::::::::::

We can specify a sample variable on which to group/organize samples along the horizontal (`x`) axis. An experimentally meaningful categorical variable is usually a good choice – in this case, the "`SampleType`" variable works much better than attempting to interpret the sample names directly (as in the previous plot):


``` r
plot_richness(GP, x="SampleType", measures=c("Chao1", "Shannon"))
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

Now suppose we wanted to use an external variable in the plot that isn’t in the `GP` dataset already – for example, a logical that indicated whether or not the samples are human-associated. First, define this new variable, `human`, as a factor (other vectors could also work; or other data you might have describing the samples).


``` r
sampleData(GP)$human <- getVariable(GP, "SampleType") %in%
  c("Feces", "Mock", "Skin", "Tongue")
```

Now tell `plot_richness` to map the new human variable on the horizontal axis, and shade the points in different color groups, according to which "`SampleType`" they belong.


``` r
plot_richness(
  GP,
  x = "human",
  color = "SampleType",
  measures = c("Chao1", "Shannon")
)
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

We can merge samples that are from the environment (`SampleType`), and make the points bigger with a ggplot2 layer. First, merge the samples.


``` r
GPst <- merge_samples(GP, "SampleType")
# repair variables that were damaged during merge (coerced to numeric)
sample_data(GPst)$SampleType <- factor(sample_names(GPst))
sample_data(GPst)$human <- as.logical(sample_data(GPst)$human)
```

Now we can plot this environment-merged version of the data. First store the default ggplot graphic as `p`, then add an additional `geom_point` layer with a large size and slight transparency.


``` r
p <- plot_richness(
  GPst,
  x = "human",
  color = "SampleType",
  measures = c("Chao1", "Shannon")
)
p1 <- p + geom_point(size = 5, alpha = 0.7)
p1
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

**More details about ggplot2**

For those interested in why this works so concisely `(p + geom_point(size=4, alpha=0.7))`, it is because the rest of the aesthetic mapping and data are contained in the ggplot object, `p`, and so is inherited in the call to the ggplot2 geometric object layer function, `geom_point`, by default since we didn’t specify alternative `aes` or `data` arguments. Although we could have if we wanted to. This perhaps sounds more confusing than it is, and I find it easier to understand by inspecting the examples I’ve shown here.

You’ll also notice that the original smaller points are still on the plot. This is because they were the first layer, and our larger points are semi-transparent. I find this kind of distracting, and doesn’t add any information or clarity. The good news is that layers can be removed from a ggplot object with standard list notation (using the dollar sign `$`).

First, check which lists are present in `p`.


``` r
p$layers
```

``` output
$geom_point
geom_point: na.rm = TRUE
stat_identity: na.rm = TRUE
position_identity 

$geom_errorbar
mapping: ymax = ~value + se, ymin = ~value - se 
geom_errorbar: na.rm = FALSE, orientation = NA, lineend = butt, width = 0.1
stat_identity: na.rm = FALSE
position_identity 
```

We can see that the first layer is the one specifying the original points, which are small. We can use negative indexing to “pop” it out, then add a new `geom_point` layer with larger point size (the following two lines).


``` r
p$layers <- p$layers[-1]
```

We can further improvde on the plots by setting the theme to something simpler.


``` r
p1 + theme_bw()
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

## Updating the data & re-formatting

Now we are going to get a bit more tedious with our data formatting and figures. Lets say on the plot above we want to actually order the `SampleType` on the legend by group the sample of human origin together and those of non-human together.

To do this - we have to first extract out the `sam_data` as its own data frame and add in a new column with the factor in this order.


``` r
# Extract the sample data as a data frame
sample_data_df <- as.data.frame(GlobalPatterns@sam_data)
# Give it a slightly different column name `SampleTypeOrigin`
sample_data_df$SampleTypeOrigin <- sample_data_df$SampleType

# Set the order for this factor
desired_order_levels <- c(
  "Feces",
  "Skin",
  "Tongue",
  "Mock",
  "Freshwater",
  "Freshwater (creek)",
  "Ocean",
  "Sediment (estuary)",
  "Soil"
)

# Now we set the factor and order
sample_data_df$SampleTypeOrigin <- factor(sample_data_df$SampleTypeOrigin, levels = desired_order_levels)
# Re-make the `GlobalPatterns` data object with this data frame and you will see we now have 9 sample variables
sample_data(GlobalPatterns) <- sample_data_df
GlobalPatterns
```

``` output
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
sample_data() Sample Data:       [ 26 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]
```

But because we remade the phyloseq object we now need to perform the same subsetting we did to generate the `GP` data object again. We also have to add in that extra `human` variable again too.


``` r
# Remove taxa with 0 reads
GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)
```

``` warning
Warning in prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns): 'prune_species' is deprecated.
Use 'prune_taxa' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` warning
Warning in speciesSums(GlobalPatterns): 'speciesSums' is deprecated.
Use 'taxa_sums' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` r
# Add in the human variable for SampleType
sampleData(GP)$human <- getVariable(GP, "SampleType") %in%
  c("Feces", "Mock", "Skin", "Tongue")
```

``` warning
Warning in getVariable(GP, "SampleType"): 'getVariable' is deprecated.
Use 'get_variable' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` warning
Warning in sampleData(`*tmp*`): 'sampleData' is deprecated.
Use 'sample_data' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` warning
Warning in `sampleData<-`(`*tmp*`, value = new("sample_data", .Data = list(: 'sampleData<-' is deprecated.
Use 'sample_data<-' instead.
See help("Deprecated") and help("phyloseq-deprecated").
```

``` r
sample_data(GP)$human <- as.logical(sample_data(GP)$human)
```

### Remaking the alpha diversity plot

Only *now* we can we remake that alpha diversity plot from before - now the change in order of the legend

``` r
p <- plot_richness(
  GP,
  x = "human",
  color = "SampleTypeOrigin",
  measures = c("Chao1", "Shannon")
)
p1 <- p + geom_point(size = 5, alpha = 0.7)
p1
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

We can then apply a manual customise colour scale for our variables, grouping by colour

``` r
p1 + scale_color_manual(values = c(
  "#9e0142",
  "#d53e4f",
  "#f46d43",
  "#252525",
  "#e6f598",
  "#abdda4",
  "#66c2a5",
  "#3288bd",
  "#5e4fa2"
) ) + theme_bw()
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

> In collaborative research, it’s almost inevitable that you’ll need to add extra variables or create new ways to subset and view your data. Different team members—each with unique backgrounds and priorities—often request specific visualizations to highlight patterns relevant to their perspective.
Rather than treating these requests as a burden, you can make them efficient and reproducible by following good practices.

When your supervisor asks for a new grouping or your collaborator wants to see results by location instead of treatment, you can adapt quickly—without rewriting your entire analysis. This approach saves time, reduces errors, and keeps your workflow transparent.


## Phylogenetic trees

The next set of figures work through some visualization pf phylogenetic trees and incorporations of sample data for annotation, examples taken from the [phyloseq tree tutorial](https://joey711.github.io/phyloseq/plot_tree-examples.html).


We want to plot trees, sometimes even bootstrap values, but notice that the node labels in the `GlobalPatterns` dataset are actually a bit strange. They look like they might be bootstrap values, but they sometimes have two decimals.


``` r
head(phy_tree(GlobalPatterns)$node.label, 10)
```

``` output
 [1] ""          "0.858.4"   "1.000.154" "0.764.3"   "0.995.2"   "1.000.2"  
 [7] "0.943.7"   "0.971.6"   "0.766"     "0.611"    
```

Could systematically remove the second decimal, but why not just take the first 4 characters instead?

``` r
phy_tree(GlobalPatterns)$node.label = substr(phy_tree(GlobalPatterns)$node.label, 1, 4)
```

Great, now that we're more happy with the node labels at least looking like bootstrap values, we can move on to using these along with other information about data mapped onto the tree graphic.

The `GlobalPatterns` dataset has many OTUs, more than we would want to try to fit on a tree graphic

``` r
ntaxa(GlobalPatterns)
```

``` output
[1] 19216
```

So, let's arbitrarily prune to just the first 50 OTUs in `GlobalPatterns`, and store this as `physeq`, which also happens to be the name for most main data parameters of function in the phyloseq package.


``` r
physeq <- prune_taxa(taxa_names(GlobalPatterns)[1:50], GlobalPatterns)
```

Now let's look at what happens with the default `plot_tree` settings.

``` r
plot_tree(physeq)
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

By default, black dots are annotated next to tips (OTUs) in the tree, one for each sample in which that OTU was observed. Some have more dots than others. Also by default, the node labels that were stored in the tree were added next to each node without any processing (although we had trimmed their length to 4 characters in the previous step).

What if we want to just see the tree with no sample points next to the tips?

``` r
plot_tree(physeq, "treeonly")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-27-1.png" style="display: block; margin: auto;" />
And what about without the node labels either?

``` r
plot_tree(physeq, "treeonly", nodeplotblank)
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-28-1.png" style="display: block; margin: auto;" />
We can adjust the way branches are rotated to make it look nicer using the `ladderize` parameter.

``` r
plot_tree(physeq, "treeonly", nodeplotblank, ladderize="left")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-29-1.png" style="display: block; margin: auto;" />

``` r
plot_tree(physeq, "treeonly", nodeplotblank, ladderize=TRUE)
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-29-2.png" style="display: block; margin: auto;" />
And what if we want to add the OTU labels next to each tip?

``` r
plot_tree(physeq, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-30-1.png" style="display: block; margin: auto;" />

Any `method` parameter argument other than `"sampledodge"` (the default) will not add dodged sample points next to the tips.

``` r
plot_tree(physeq, "anythingelse")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-31-1.png" style="display: block; margin: auto;" />

### Mapping Variables in Data

In the default argument to `method`, `"sampledodge"`, a point is added next to each OTU tip in the tree for every sample in which that OTU was observed. We can then map certain aesthetic features of these points to variables in our data.

#### Color

Color is one of the most useful aesthetics in tree graphics when they are complicated. Color can be mapped to either taxonomic ranks or sample covariates. For instance, we can map color to the type of sample collected (environmental location).

``` r
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left", color="SampleType")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

**Keeping colours consistant!**

Now that you are intermediate R users you would be picking up on the different colour scales used to map against the same variable. It is best practice to keep these consistent across the different types of figures. So let's regenerate the figure mapping the `SampleTypeOrigin` variable to the same set of colours we used above for the alpha diversity plot


``` r
plot_tree(physeq,
  nodelabf = nodeplotboot(),
  ladderize = "left",
  color = "SampleTypeOrigin") + scale_color_manual(
    values = c(
      "#9e0142",
      "#d53e4f",
      "#f46d43",
      "#252525",
      "#e6f598",
      "#abdda4",
      "#66c2a5",
      "#3288bd",
      "#5e4fa2"
    )
  ) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

You can make this even more efficent by including the colouring of a variable in your global R chunk for the set up and then refer to it using the value - in this case `SampleTypeOriginCols`


``` r
SampleTypeOriginCols <- c(
      "#9e0142",
      "#d53e4f",
      "#f46d43",
      "#252525",
      "#e6f598",
      "#abdda4",
      "#66c2a5",
      "#3288bd",
      "#5e4fa2"
    )
plot_tree(physeq,
  nodelabf = nodeplotboot(),
  ladderize = "left",
  color = "SampleTypeOrigin") + scale_color_manual(
    values = SampleTypeOriginCols)
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-34-1.png" style="display: block; margin: auto;" />


We can also map color to taxonomic class.

``` r
plot_tree(physeq, nodelabf=nodeplotboot(), ladderize="left", color="Class")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-35-1.png" style="display: block; margin: auto;" />
This time we don't mind the differrent colours because we are mapping it against taxonomy

#### Shape

You can also map a variable to point shape if it has 6 or fewer categories, and this can be done even when color is also mapped. Here we map shape to taxonomic class so that we can still indicate it in the graphic while also mapping `SampleType` to point color.

``` r
plot_tree(
  physeq,
  nodelabf = nodeplotboot(),
  ladderize = "left",
  color = "SampleTypeOrigin",
  shape = "Class"
) + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-36-1.png" style="display: block; margin: auto;" />

### Node labels

One of the most common reasons to label nodes is to add confidence measures, often a bootstrap value, to the nodes of the tree. The following graphics show different ways of doing this (labels are added by default if present in your tree).

``` r
# The default
plot_tree(physeq, color = "SampleTypeOrigin", ladderize = "left") + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-37-1.png" style="display: block; margin: auto;" />

``` r
# Special bootstrap label
plot_tree(physeq,
  nodelabf = nodeplotboot(),
  color = "SampleTypeOrigin",
  ladderize = "left") + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-37-2.png" style="display: block; margin: auto;" />

``` r
# Special bootstrap label with alternative thresholds
plot_tree(
  physeq,
  nodelabf = nodeplotboot(80, 0, 3),
  color = "SampleTypeOrigin",
  ladderize = "left"
) + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-37-3.png" style="display: block; margin: auto;" />

### Tip labels

- **label.tips** - The `label.tips` parameter controls labeling of tree tips (AKA leaves).  Default is NULL, indicating that no tip labels will be printed. If `"taxa_names"` is a special argument resulting in the OTU name (try `taxa_names` function) being labelled next to the leaves or next to the set of points that label the leaves. Alternatively, if your data object contains a `tax_table`, then one of the rank names (from `rank_names(physeq)`) can be provided, and the classification of each OTU at that rank will be labeled instead.
- **text.size** - A positive numeric argument indicating the ggplot2 size parameter for the taxa labels. Default is `NULL`. If left as `NULL`, this function will automatically calculate a (hopefully) optimal text size given the size constraints posed by the tree itself (for a vertical tree). This argument is included mainly in case the automatically-calculated size is wrong and you want to change it. Note that this parameter is only meaningful if `label.tips` is not `NULL`


``` r
plot_tree(
  physeq,
  nodelabf = nodeplotboot(80, 0, 3),
  color = "SampleTypeOrigin",
  label.tips = "taxa_names",
  ladderize = "left"
) + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-38-1.png" style="display: block; margin: auto;" />

### Radial Tree

Making a radial tree is easy with ggplot2, simply recognizing that our vertically-oriented tree is a cartesian mapping of the data to a graphic -- and that a radial tree is the same mapping, but with polar coordinates instead.


``` r
data(esophagus)
plot_tree(esophagus, color = "Sample", ladderize = "left") + coord_polar(theta =
    "y")
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-39-1.png" style="display: block; margin: auto;" />

The `GlobalPatterns` dataset has additional data we can map, so we will re-do some preliminary data loading/trimming to make this radial-tree example self contained, and then show the same plot as above.


``` r
plot_tree(
  physeq,
  nodelabf = nodeplotboot(60, 60, 3),
  color = "SampleTypeOrigin",
  shape = "Class",
  ladderize = "left"
) + coord_polar(theta = "y") + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-40-1.png" style="display: block; margin: auto;" />


:::::::::::::::::::::::::::::::::::::::  challenge

### More Examples with the Global Patterns dataset

Try exploring more options for generating trees but subsetting specific taxa. Starting the Archaea (*hint: Archaea is a kindgom*).

:::::::::::::::  solution

### Solution


``` r
gpa <- subset_taxa(GlobalPatterns, Kingdom=="Archaea")
```

Identify if the number of different Archaeal taxa is suitable for a tree


``` r
ntaxa(gpa)
```

``` output
[1] 208
```

Visualise tree

``` r
plot_tree(
  gpa,
  nodelabf = nodeplotboot(80, 0, 3),
  color = "SampleTypeOrigin",
  label.tips = "taxa_names",
  ladderize = "left"
) + scale_color_manual(values = SampleTypeOriginCols) 
```

<img src="fig/06-application-microbiome-rendered-unnamed-chunk-43-1.png" style="display: block; margin: auto;" />

:::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::::


https://microbiome.github.io/OMA/docs/devel/

:::::::::::::::::::::::::::::::::::::::::  callout

### Want more bioconductor and microbiome (and related) analysis?

- [Orchestrating Microbiome Analysis with Bioconductor](https://microbiome.github.io/OMA/docs/devel/) and related packages from the [miaverse](https://microbiome.github.io/)
- [*microbiome* package tutorials](https://microbiome.github.io/tutorials/)
- [*phyloseq* package tutorials](https://joey711.github.io/phyloseq/)

::::::::::::::::::::::::::::::::::::::::::::::::::

---

:::::::::::::::::::::::::::::::::::::::: keypoints

- Bespoke often domain specific data classes are often built to serve a specific set of data, however by utilising fundamental R principals such as S4 class structures transferring formats between packages for analysis is much easier.
- Look for widely accepted data classes in your field before formatting your data for a specific package early on.
- Packages that generate [*ggplot2*](https://ggplot2.tidyverse.org/) compatible figures allow for *easier* manipulations and customise using familar syntax
- Get feedback early on from preliminary analysis, this will help you improve your data and reporting for generating downstream analysis and figures.
- Get comfortabe with the data lifecycle and **re**-producing analysis and figures!

::::::::::::::::::::::::::::::::::::::::::::::::::
