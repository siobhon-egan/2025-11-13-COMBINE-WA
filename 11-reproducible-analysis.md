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
install.packages("DT")

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


``` warning
Warning: `includeHTML()` was provided a `path` that appears to be a complete HTML document.
✖ Path: files/microbiome-data-report.html
ℹ Use `tags$iframe()` to include an HTML document. You can either ensure `path` is accessible in your app or document (see e.g. `shiny::addResourcePath()`) and pass the relative path to the `src` argument. Or you can read the contents of `path` and pass the contents to `srcdoc`.
```

<!--html_preserve--><!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en"><head>

<meta charset="utf-8">
<meta name="generator" content="quarto-1.7.23">

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes">

<meta name="author" content="Siobhon Egan">
<meta name="dcterms.date" content="2025-11-13">

<title>Microbiome data report for Supervisor meeting</title>
<style>
code{white-space: pre-wrap;}
span.smallcaps{font-variant: small-caps;}
div.columns{display: flex; gap: min(4vw, 1.5em);}
div.column{flex: auto; overflow-x: auto;}
div.hanging-indent{margin-left: 1.5em; text-indent: -1.5em;}
ul.task-list{list-style: none;}
ul.task-list li input[type="checkbox"] {
  width: 0.8em;
  margin: 0 0.8em 0.2em -1em; /* quarto-specific, see https://github.com/quarto-dev/quarto-cli/issues/4556 */ 
  vertical-align: middle;
}
/* CSS for syntax highlighting */
html { -webkit-text-size-adjust: 100%; }
pre > code.sourceCode { white-space: pre; position: relative; }
pre > code.sourceCode > span { display: inline-block; line-height: 1.25; }
pre > code.sourceCode > span:empty { height: 1.2em; }
.sourceCode { overflow: visible; }
code.sourceCode > span { color: inherit; text-decoration: inherit; }
div.sourceCode { margin: 1em 0; }
pre.sourceCode { margin: 0; }
@media screen {
div.sourceCode { overflow: auto; }
}
@media print {
pre > code.sourceCode { white-space: pre-wrap; }
pre > code.sourceCode > span { text-indent: -5em; padding-left: 5em; }
}
pre.numberSource code
  { counter-reset: source-line 0; }
pre.numberSource code > span
  { position: relative; left: -4em; counter-increment: source-line; }
pre.numberSource code > span > a:first-child::before
  { content: counter(source-line);
    position: relative; left: -1em; text-align: right; vertical-align: baseline;
    border: none; display: inline-block;
    -webkit-touch-callout: none; -webkit-user-select: none;
    -khtml-user-select: none; -moz-user-select: none;
    -ms-user-select: none; user-select: none;
    padding: 0 4px; width: 4em;
  }
pre.numberSource { margin-left: 3em;  padding-left: 4px; }
div.sourceCode
  {   }
@media screen {
pre > code.sourceCode > span > a:first-child::before { text-decoration: underline; }
}
</style>


<script src="microbiome-data-report_files/libs/clipboard/clipboard.min.js"></script>
<script src="microbiome-data-report_files/libs/quarto-html/quarto.js" type="module"></script>
<script src="microbiome-data-report_files/libs/quarto-html/tabsets/tabsets.js" type="module"></script>
<script src="microbiome-data-report_files/libs/quarto-html/popper.min.js"></script>
<script src="microbiome-data-report_files/libs/quarto-html/tippy.umd.min.js"></script>
<script src="microbiome-data-report_files/libs/quarto-html/anchor.min.js"></script>
<link href="microbiome-data-report_files/libs/quarto-html/tippy.css" rel="stylesheet">
<link href="microbiome-data-report_files/libs/quarto-html/quarto-syntax-highlighting-226bd0f977fa82dfae4534cac220d79a.css" rel="stylesheet" id="quarto-text-highlighting-styles">
<script src="microbiome-data-report_files/libs/bootstrap/bootstrap.min.js"></script>
<link href="microbiome-data-report_files/libs/bootstrap/bootstrap-icons.css" rel="stylesheet">
<link href="microbiome-data-report_files/libs/bootstrap/bootstrap-be9206a58a27969b3ea3001da8f26fef.min.css" rel="stylesheet" append-hash="true" id="quarto-bootstrap" data-mode="light">


</head>

<body class="fullcontent quarto-light">

<div id="quarto-content" class="page-columns page-rows-contents page-layout-article">

<main class="content" id="quarto-document-content">

<header id="title-block-header" class="quarto-title-block default">
<div class="quarto-title">
<h1 class="title">Microbiome data report for Supervisor meeting</h1>
</div>



<div class="quarto-title-meta">

    <div>
    <div class="quarto-title-meta-heading">Author</div>
    <div class="quarto-title-meta-contents">
             <p>Siobhon Egan </p>
          </div>
  </div>
    
    <div>
    <div class="quarto-title-meta-heading">Published</div>
    <div class="quarto-title-meta-contents">
      <p class="date">November 13, 2025</p>
    </div>
  </div>
  
    
  </div>
  


</header>


<p>The following datasets and analysis have been used from tutorials available in the <a href="">microbiome package</a>, for <a href="https://microbiome.github.io/tutorials/Heatmap.html">heatmaps</a></p>
<section id="information" class="level2">
<h2 class="anchored" data-anchor-id="information">Information</h2>
<p>Prepare a summary table of the participants in your dataset</p>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb1"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb1-1"><a href="#cb1-1" aria-hidden="true" tabindex="-1"></a>pseq<span class="sc">@</span>sam_data <span class="sc">%&gt;%</span></span>
<span id="cb1-2"><a href="#cb1-2" aria-hidden="true" tabindex="-1"></a>  <span class="fu">group_by</span>(group, sex, time) <span class="sc">%&gt;%</span></span>
<span id="cb1-3"><a href="#cb1-3" aria-hidden="true" tabindex="-1"></a>  <span class="fu">summarise</span>(<span class="at">n =</span> <span class="fu">n</span>())</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output cell-output-stdout">
<pre><code># A tibble: 8 × 4
# Groups:   group, sex [4]
  group   sex     time     n
  &lt;fct&gt;   &lt;fct&gt;  &lt;int&gt; &lt;int&gt;
1 LGG     female     1     5
2 LGG     female     2     5
3 LGG     male       1     3
4 LGG     male       2     3
5 Placebo female     1    10
6 Placebo female     2    10
7 Placebo male       1     4
8 Placebo male       2     4</code></pre>
</div>
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb3"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb3-1"><a href="#cb3-1" aria-hidden="true" tabindex="-1"></a>pseq<span class="sc">@</span>sam_data <span class="sc">%&gt;%</span></span>
<span id="cb3-2"><a href="#cb3-2" aria-hidden="true" tabindex="-1"></a>  <span class="fu">group_by</span>(group, sex) <span class="sc">%&gt;%</span></span>
<span id="cb3-3"><a href="#cb3-3" aria-hidden="true" tabindex="-1"></a>  <span class="fu">summarise</span>(<span class="at">n =</span> <span class="fu">n</span>())</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output cell-output-stdout">
<pre><code># A tibble: 4 × 3
# Groups:   group [2]
  group   sex        n
  &lt;fct&gt;   &lt;fct&gt;  &lt;int&gt;
1 LGG     female    10
2 LGG     male       6
3 Placebo female    20
4 Placebo male       8</code></pre>
</div>
</div>
</section>
<section id="alpha-diversity" class="level2">
<h2 class="anchored" data-anchor-id="alpha-diversity">Alpha diversity</h2>
</section>
<section id="beta-diversity" class="level2">
<h2 class="anchored" data-anchor-id="beta-diversity">Beta diversity</h2>
</section>
<section id="heatmap-of-abundance" class="level2">
<h2 class="anchored" data-anchor-id="heatmap-of-abundance">Heatmap of abundance</h2>
<p>Abundance of microbial taxa (using the Z transformed data)</p>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb5"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb5-1"><a href="#cb5-1" aria-hidden="true" tabindex="-1"></a><span class="co"># Plot the abundances heatmap, round values to 2 decimals</span></span>
<span id="cb5-2"><a href="#cb5-2" aria-hidden="true" tabindex="-1"></a>dfm <span class="ot">&lt;-</span> <span class="fu">melt</span>(<span class="fu">round</span>(<span class="fu">abundances</span>(pseqz), <span class="dv">1</span>))</span>
<span id="cb5-3"><a href="#cb5-3" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb5-4"><a href="#cb5-4" aria-hidden="true" tabindex="-1"></a><span class="fu">colnames</span>(dfm) <span class="ot">&lt;-</span> <span class="fu">c</span>(<span class="st">"Taxa"</span>, <span class="st">"Sample"</span>, <span class="st">"value"</span>)</span>
<span id="cb5-5"><a href="#cb5-5" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb5-6"><a href="#cb5-6" aria-hidden="true" tabindex="-1"></a><span class="fu">heat</span>(dfm, <span class="st">"Taxa"</span>, <span class="st">"Sample"</span>, <span class="st">"value"</span>) <span class="sc">+</span></span>
<span id="cb5-7"><a href="#cb5-7" aria-hidden="true" tabindex="-1"></a>  <span class="fu">theme</span>(<span class="at">text=</span><span class="fu">element_text</span>(<span class="at">size=</span><span class="dv">10</span>), </span>
<span id="cb5-8"><a href="#cb5-8" aria-hidden="true" tabindex="-1"></a>        <span class="at">axis.text.x =</span> <span class="fu">element_text</span>(<span class="at">angle =</span> <span class="dv">90</span>, <span class="at">hjust =</span> <span class="dv">1</span>),</span>
<span id="cb5-9"><a href="#cb5-9" aria-hidden="true" tabindex="-1"></a>        <span class="at">legend.key.size =</span> <span class="fu">unit</span>(<span class="fl">1.2</span>, <span class="st">"cm"</span>))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="microbiome-data-report_files/figure-html/unnamed-chunk-3-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
</section>
<section id="correlation-between-microbes-and-lipids" class="level2">
<h2 class="anchored" data-anchor-id="correlation-between-microbes-and-lipids">Correlation between microbes and lipids</h2>
<p>Cross-correlate columns of two data sets from related to microbiome and blood serum lipids associations (<a href="https://peerj.com/articles/32/">PeerJ 1:e32</a>).</p>
<p>The function returns correlations, raw p-values, and fdr estimates (not strictly proper as the comparisons are not independent). Keep only those elements that have at least only one significant correlation (n.signif):</p>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb6"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb6-1"><a href="#cb6-1" aria-hidden="true" tabindex="-1"></a><span class="co"># Load example data </span></span>
<span id="cb6-2"><a href="#cb6-2" aria-hidden="true" tabindex="-1"></a>otu <span class="ot">&lt;-</span> peerj32<span class="sc">$</span>microbes </span>
<span id="cb6-3"><a href="#cb6-3" aria-hidden="true" tabindex="-1"></a>lipids <span class="ot">&lt;-</span> peerj32<span class="sc">$</span>lipids </span>
<span id="cb6-4"><a href="#cb6-4" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-5"><a href="#cb6-5" aria-hidden="true" tabindex="-1"></a><span class="co"># Define data sets to cross-correlate</span></span>
<span id="cb6-6"><a href="#cb6-6" aria-hidden="true" tabindex="-1"></a>x <span class="ot">&lt;-</span> <span class="fu">log10</span>(otu) <span class="co"># OTU Log10 (44 samples x 130 genera)</span></span>
<span id="cb6-7"><a href="#cb6-7" aria-hidden="true" tabindex="-1"></a>y <span class="ot">&lt;-</span> <span class="fu">as.matrix</span>(lipids) <span class="co"># Lipids (44 samples x 389 lipids)</span></span>
<span id="cb6-8"><a href="#cb6-8" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-9"><a href="#cb6-9" aria-hidden="true" tabindex="-1"></a><span class="co"># Cross correlate data sets</span></span>
<span id="cb6-10"><a href="#cb6-10" aria-hidden="true" tabindex="-1"></a>correlations <span class="ot">&lt;-</span> <span class="fu">associate</span>(x, y, <span class="at">method =</span> <span class="st">"spearman"</span>, <span class="at">mode =</span> <span class="st">"matrix"</span>, <span class="at">p.adj.threshold =</span> <span class="fl">0.05</span>, <span class="at">n.signif =</span> <span class="dv">1</span>)</span>
<span id="cb6-11"><a href="#cb6-11" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-12"><a href="#cb6-12" aria-hidden="true" tabindex="-1"></a><span class="co"># Or, alternatively, the same output is also available in a handy table format</span></span>
<span id="cb6-13"><a href="#cb6-13" aria-hidden="true" tabindex="-1"></a>correlation.table <span class="ot">&lt;-</span> <span class="fu">associate</span>(x, y, <span class="at">method =</span> <span class="st">"spearman"</span>, <span class="at">mode =</span> <span class="st">"table"</span>, <span class="at">p.adj.threshold =</span> <span class="fl">0.05</span>, <span class="at">n.signif =</span> <span class="dv">1</span>)</span>
<span id="cb6-14"><a href="#cb6-14" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb6-15"><a href="#cb6-15" aria-hidden="true" tabindex="-1"></a><span class="fu">kable</span>(<span class="fu">head</span>(correlation.table))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output-display">
<table class="caption-top table table-sm table-striped small">
<colgroup>
<col style="width: 7%">
<col style="width: 45%">
<col style="width: 15%">
<col style="width: 17%">
<col style="width: 14%">
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: left;">X1</th>
<th style="text-align: left;">X2</th>
<th style="text-align: right;">Correlation</th>
<th style="text-align: right;">p.adj</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">1648</td>
<td style="text-align: left;">Ruminococcus gnavus et rel.</td>
<td style="text-align: left;">TG(54:5).2</td>
<td style="text-align: right;">0.7164958</td>
<td style="text-align: right;">0.0022842</td>
</tr>
<tr class="even">
<td style="text-align: left;">384</td>
<td style="text-align: left;">Moraxellaceae</td>
<td style="text-align: left;">PC(40:3e)</td>
<td style="text-align: right;">-0.6941863</td>
<td style="text-align: right;">0.0029225</td>
</tr>
<tr class="odd">
<td style="text-align: left;">1829</td>
<td style="text-align: left;">Uncultured Bacteroidetes</td>
<td style="text-align: left;">TG(56:2).1</td>
<td style="text-align: right;">-0.6987375</td>
<td style="text-align: right;">0.0029225</td>
</tr>
<tr class="even">
<td style="text-align: left;">349</td>
<td style="text-align: left;">Lactobacillus plantarum et rel.</td>
<td style="text-align: left;">PC(40:3)</td>
<td style="text-align: right;">-0.6877976</td>
<td style="text-align: right;">0.0031520</td>
</tr>
<tr class="odd">
<td style="text-align: left;">1198</td>
<td style="text-align: left;">Ruminococcus gnavus et rel.</td>
<td style="text-align: left;">TG(52:5)</td>
<td style="text-align: right;">0.6806216</td>
<td style="text-align: right;">0.0037518</td>
</tr>
<tr class="even">
<td style="text-align: left;">264</td>
<td style="text-align: left;">Moraxellaceae</td>
<td style="text-align: left;">PC(38:4).1</td>
<td style="text-align: right;">-0.6700504</td>
<td style="text-align: right;">0.0038414</td>
</tr>
</tbody>
</table>
</div>
</div>
<section id="association-heatmaps" class="level3">
<h3 class="anchored" data-anchor-id="association-heatmaps">Association heatmaps</h3>
<p>Rearrange the data and plot the heatmap and mark significant correlations with stars to reproduce microbiota-lipidome heatmap from <a href="https://peerj.com/articles/32/">Lahti et al.&nbsp;PeerJ (2013)</a> (the ordering of rows and columns may be different):</p>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb7"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb7-1"><a href="#cb7-1" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> <span class="fu">heat</span>(correlation.table, <span class="st">"X1"</span>, <span class="st">"X2"</span>, </span>
<span id="cb7-2"><a href="#cb7-2" aria-hidden="true" tabindex="-1"></a>          <span class="at">fill =</span> <span class="st">"Correlation"</span>, </span>
<span id="cb7-3"><a href="#cb7-3" aria-hidden="true" tabindex="-1"></a>          <span class="at">star =</span> <span class="st">"p.adj"</span>, </span>
<span id="cb7-4"><a href="#cb7-4" aria-hidden="true" tabindex="-1"></a>          <span class="at">p.adj.threshold =</span> <span class="fl">0.05</span>) </span>
<span id="cb7-5"><a href="#cb7-5" aria-hidden="true" tabindex="-1"></a>p <span class="sc">+</span> <span class="fu">theme</span>(<span class="at">text=</span><span class="fu">element_text</span>(<span class="at">size=</span><span class="dv">10</span>), </span>
<span id="cb7-6"><a href="#cb7-6" aria-hidden="true" tabindex="-1"></a>          <span class="at">axis.text.x =</span> <span class="fu">element_text</span>(<span class="at">angle =</span> <span class="dv">90</span>, <span class="at">hjust =</span> <span class="dv">1</span>), </span>
<span id="cb7-7"><a href="#cb7-7" aria-hidden="true" tabindex="-1"></a>          <span class="at">legend.key.size =</span> <span class="fu">unit</span>(<span class="fl">1.3</span>, <span class="st">"cm"</span>))</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="microbiome-data-report_files/figure-html/unnamed-chunk-5-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
</section>
<section id="heatmap-of-significant-correlations" class="level3">
<h3 class="anchored" data-anchor-id="heatmap-of-significant-correlations">Heatmap of significant correlations</h3>
<p>This time we are modifying the plot to only include significant correlations with a threshold of P &lt; 0.05. The cells marked with a white cross indicate most significant correlations with a significance of P &lt; 0.02.</p>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb8"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb8-1"><a href="#cb8-1" aria-hidden="true" tabindex="-1"></a><span class="co"># Order the rows and columns with levels argument if needed:</span></span>
<span id="cb8-2"><a href="#cb8-2" aria-hidden="true" tabindex="-1"></a>correlation.table<span class="sc">$</span>X1 <span class="ot">&lt;-</span> <span class="fu">factor</span>(correlation.table<span class="sc">$</span>X1, <span class="at">levels =</span> <span class="fu">unique</span>(<span class="fu">as.character</span>(correlation.table<span class="sc">$</span>X1)))</span>
<span id="cb8-3"><a href="#cb8-3" aria-hidden="true" tabindex="-1"></a>correlation.table<span class="sc">$</span>X2 <span class="ot">&lt;-</span> <span class="fu">factor</span>(correlation.table<span class="sc">$</span>X2, <span class="at">levels =</span> <span class="fu">unique</span>(<span class="fu">as.character</span>(correlation.table<span class="sc">$</span>X2)))</span>
<span id="cb8-4"><a href="#cb8-4" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-5"><a href="#cb8-5" aria-hidden="true" tabindex="-1"></a><span class="co"># Set black-and-white theme</span></span>
<span id="cb8-6"><a href="#cb8-6" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(ggplot2)</span>
<span id="cb8-7"><a href="#cb8-7" aria-hidden="true" tabindex="-1"></a><span class="fu">theme_set</span>(<span class="fu">theme_bw</span>())</span>
<span id="cb8-8"><a href="#cb8-8" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-9"><a href="#cb8-9" aria-hidden="true" tabindex="-1"></a><span class="co"># Pick only the correlations with p.adj&lt;0.05</span></span>
<span id="cb8-10"><a href="#cb8-10" aria-hidden="true" tabindex="-1"></a><span class="co"># Note: this will leave other cells empty</span></span>
<span id="cb8-11"><a href="#cb8-11" aria-hidden="true" tabindex="-1"></a><span class="fu">library</span>(dplyr)</span>
<span id="cb8-12"><a href="#cb8-12" aria-hidden="true" tabindex="-1"></a>subtable <span class="ot">&lt;-</span> <span class="fu">filter</span>(correlation.table, p.adj <span class="sc">&lt;</span> <span class="fl">0.05</span>)</span>
<span id="cb8-13"><a href="#cb8-13" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-14"><a href="#cb8-14" aria-hidden="true" tabindex="-1"></a><span class="co"># Arrange the figure</span></span>
<span id="cb8-15"><a href="#cb8-15" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> <span class="fu">ggplot</span>(subtable, <span class="fu">aes</span>(<span class="at">x =</span> X1, <span class="at">y =</span> X2, <span class="at">fill =</span> Correlation))</span>
<span id="cb8-16"><a href="#cb8-16" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> p <span class="sc">+</span> <span class="fu">geom_tile</span>() </span>
<span id="cb8-17"><a href="#cb8-17" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> p <span class="sc">+</span> <span class="fu">scale_fill_gradientn</span>(<span class="st">"Correlation"</span>, </span>
<span id="cb8-18"><a href="#cb8-18" aria-hidden="true" tabindex="-1"></a>                       <span class="at">breaks =</span> <span class="fu">seq</span>(<span class="at">from =</span> <span class="sc">-</span><span class="dv">1</span>, <span class="at">to =</span> <span class="dv">1</span>, <span class="at">by =</span> <span class="fl">0.2</span>), </span>
<span id="cb8-19"><a href="#cb8-19" aria-hidden="true" tabindex="-1"></a>                   <span class="at">colours =</span> <span class="fu">c</span>(<span class="st">"darkblue"</span>, <span class="st">"blue"</span>, <span class="st">"white"</span>, <span class="st">"red"</span>, <span class="st">"darkred"</span>), </span>
<span id="cb8-20"><a href="#cb8-20" aria-hidden="true" tabindex="-1"></a>                   <span class="at">limits =</span> <span class="fu">c</span>(<span class="sc">-</span><span class="dv">1</span>,<span class="dv">1</span>)) </span>
<span id="cb8-21"><a href="#cb8-21" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-22"><a href="#cb8-22" aria-hidden="true" tabindex="-1"></a><span class="co"># Polish texts</span></span>
<span id="cb8-23"><a href="#cb8-23" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> p <span class="sc">+</span> <span class="fu">theme</span>(<span class="at">axis.text.x=</span><span class="fu">element_text</span>(<span class="at">angle =</span> <span class="dv">90</span>, <span class="at">hjust=</span><span class="dv">1</span>, <span class="at">face =</span> <span class="st">"italic"</span>),</span>
<span id="cb8-24"><a href="#cb8-24" aria-hidden="true" tabindex="-1"></a>               <span class="at">axis.text.y=</span><span class="fu">element_text</span>(<span class="at">size =</span> <span class="dv">8</span>))</span>
<span id="cb8-25"><a href="#cb8-25" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> p <span class="sc">+</span> <span class="fu">xlab</span>(<span class="st">""</span>) <span class="sc">+</span> <span class="fu">ylab</span>(<span class="st">""</span>)</span>
<span id="cb8-26"><a href="#cb8-26" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-27"><a href="#cb8-27" aria-hidden="true" tabindex="-1"></a><span class="co"># Mark the most significant cells with stars using threshold of P &lt; 0.02</span></span>
<span id="cb8-28"><a href="#cb8-28" aria-hidden="true" tabindex="-1"></a>p <span class="ot">&lt;-</span> p <span class="sc">+</span> <span class="fu">geom_text</span>(<span class="at">data =</span> <span class="fu">subset</span>(correlation.table, p.adj <span class="sc">&lt;</span> <span class="fl">0.02</span>), </span>
<span id="cb8-29"><a href="#cb8-29" aria-hidden="true" tabindex="-1"></a>               <span class="fu">aes</span>(<span class="at">x =</span> X1, <span class="at">y =</span> X2, <span class="at">label =</span> <span class="st">"+"</span>), <span class="at">col =</span> <span class="st">"white"</span>, <span class="at">size =</span> <span class="dv">5</span>)</span>
<span id="cb8-30"><a href="#cb8-30" aria-hidden="true" tabindex="-1"></a></span>
<span id="cb8-31"><a href="#cb8-31" aria-hidden="true" tabindex="-1"></a><span class="co"># Plot</span></span>
<span id="cb8-32"><a href="#cb8-32" aria-hidden="true" tabindex="-1"></a><span class="fu">print</span>(p)</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output-display">
<div>
<figure class="figure">
<p><img src="microbiome-data-report_files/figure-html/unnamed-chunk-6-1.png" class="img-fluid figure-img" width="672"></p>
</figure>
</div>
</div>
</div>
<div class="cell">
<details class="code-fold">
<summary>Show the code</summary>
<div class="sourceCode cell-code" id="cb9"><pre class="sourceCode r code-with-copy"><code class="sourceCode r"><span id="cb9-1"><a href="#cb9-1" aria-hidden="true" tabindex="-1"></a><span class="fu">sessionInfo</span>()</span></code><button title="Copy to Clipboard" class="code-copy-button"><i class="bi"></i></button></pre></div>
</details>
<div class="cell-output cell-output-stdout">
<pre><code>R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Sonoma 14.7

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: Australia/Perth
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] gtsummary_2.4.0   knitr_1.50        reshape2_1.4.4    dplyr_1.1.4      
[5] microbiome_1.30.0 ggplot2_4.0.0     phyloseq_1.52.0  

loaded via a namespace (and not attached):
 [1] gtable_0.3.6            xfun_0.53               htmlwidgets_1.6.4      
 [4] rhdf5_2.52.1            Biobase_2.68.0          lattice_0.22-7         
 [7] rhdf5filters_1.20.0     vctrs_0.6.5             tools_4.5.1            
[10] generics_0.1.4          biomformat_1.36.0       stats4_4.5.1           
[13] parallel_4.5.1          tibble_3.3.0            cluster_2.1.8.1        
[16] pkgconfig_2.0.3         Matrix_1.7-3            data.table_1.17.8      
[19] RColorBrewer_1.1-3      S7_0.2.0                S4Vectors_0.46.0       
[22] lifecycle_1.0.4         GenomeInfoDbData_1.2.14 compiler_4.5.1         
[25] farver_2.1.2            stringr_1.5.2           Biostrings_2.76.0      
[28] codetools_0.2-20        permute_0.9-8           GenomeInfoDb_1.44.3    
[31] htmltools_0.5.8.1       yaml_2.3.10             tidyr_1.3.1            
[34] pillar_1.11.1           crayon_1.5.3            MASS_7.3-65            
[37] vegan_2.7-2             iterators_1.0.14        foreach_1.5.2          
[40] nlme_3.1-168            tidyselect_1.2.1        digest_0.6.37          
[43] Rtsne_0.17              stringi_1.8.7           purrr_1.1.0            
[46] splines_4.5.1           ade4_1.7-23             fastmap_1.2.0          
[49] grid_4.5.1              cli_3.6.5               magrittr_2.0.4         
[52] utf8_1.2.6              survival_3.8-3          ape_5.8-1              
[55] withr_3.0.2             UCSC.utils_1.4.0        scales_1.4.0           
[58] rmarkdown_2.30          XVector_0.48.0          httr_1.4.7             
[61] igraph_2.2.0            multtest_2.64.0         evaluate_1.0.5         
[64] IRanges_2.42.0          mgcv_1.9-3              rlang_1.1.6            
[67] Rcpp_1.1.0              glue_1.8.0              BiocGenerics_0.54.1    
[70] rstudioapi_0.17.1       jsonlite_2.0.0          R6_2.6.1               
[73] Rhdf5lib_1.30.0         plyr_1.8.9             </code></pre>
</div>
</div>
</section>
</section>

</main>
<!-- /main column -->
<script id="quarto-html-after-body" type="application/javascript">
  window.document.addEventListener("DOMContentLoaded", function (event) {
    const icon = "";
    const anchorJS = new window.AnchorJS();
    anchorJS.options = {
      placement: 'right',
      icon: icon
    };
    anchorJS.add('.anchored');
    const isCodeAnnotation = (el) => {
      for (const clz of el.classList) {
        if (clz.startsWith('code-annotation-')) {                     
          return true;
        }
      }
      return false;
    }
    const onCopySuccess = function(e) {
      // button target
      const button = e.trigger;
      // don't keep focus
      button.blur();
      // flash "checked"
      button.classList.add('code-copy-button-checked');
      var currentTitle = button.getAttribute("title");
      button.setAttribute("title", "Copied!");
      let tooltip;
      if (window.bootstrap) {
        button.setAttribute("data-bs-toggle", "tooltip");
        button.setAttribute("data-bs-placement", "left");
        button.setAttribute("data-bs-title", "Copied!");
        tooltip = new bootstrap.Tooltip(button, 
          { trigger: "manual", 
            customClass: "code-copy-button-tooltip",
            offset: [0, -8]});
        tooltip.show();    
      }
      setTimeout(function() {
        if (tooltip) {
          tooltip.hide();
          button.removeAttribute("data-bs-title");
          button.removeAttribute("data-bs-toggle");
          button.removeAttribute("data-bs-placement");
        }
        button.setAttribute("title", currentTitle);
        button.classList.remove('code-copy-button-checked');
      }, 1000);
      // clear code selection
      e.clearSelection();
    }
    const getTextToCopy = function(trigger) {
        const codeEl = trigger.previousElementSibling.cloneNode(true);
        for (const childEl of codeEl.children) {
          if (isCodeAnnotation(childEl)) {
            childEl.remove();
          }
        }
        return codeEl.innerText;
    }
    const clipboard = new window.ClipboardJS('.code-copy-button:not([data-in-quarto-modal])', {
      text: getTextToCopy
    });
    clipboard.on('success', onCopySuccess);
    if (window.document.getElementById('quarto-embedded-source-code-modal')) {
      const clipboardModal = new window.ClipboardJS('.code-copy-button[data-in-quarto-modal]', {
        text: getTextToCopy,
        container: window.document.getElementById('quarto-embedded-source-code-modal')
      });
      clipboardModal.on('success', onCopySuccess);
    }
      var localhostRegex = new RegExp(/^(?:http|https):\/\/localhost\:?[0-9]*\//);
      var mailtoRegex = new RegExp(/^mailto:/);
        var filterRegex = new RegExp('/' + window.location.host + '/');
      var isInternal = (href) => {
          return filterRegex.test(href) || localhostRegex.test(href) || mailtoRegex.test(href);
      }
      // Inspect non-navigation links and adorn them if external
     var links = window.document.querySelectorAll('a[href]:not(.nav-link):not(.navbar-brand):not(.toc-action):not(.sidebar-link):not(.sidebar-item-toggle):not(.pagination-link):not(.no-external):not([aria-hidden]):not(.dropdown-item):not(.quarto-navigation-tool):not(.about-link)');
      for (var i=0; i<links.length; i++) {
        const link = links[i];
        if (!isInternal(link.href)) {
          // undo the damage that might have been done by quarto-nav.js in the case of
          // links that we want to consider external
          if (link.dataset.originalHref !== undefined) {
            link.href = link.dataset.originalHref;
          }
        }
      }
    function tippyHover(el, contentFn, onTriggerFn, onUntriggerFn) {
      const config = {
        allowHTML: true,
        maxWidth: 500,
        delay: 100,
        arrow: false,
        appendTo: function(el) {
            return el.parentElement;
        },
        interactive: true,
        interactiveBorder: 10,
        theme: 'quarto',
        placement: 'bottom-start',
      };
      if (contentFn) {
        config.content = contentFn;
      }
      if (onTriggerFn) {
        config.onTrigger = onTriggerFn;
      }
      if (onUntriggerFn) {
        config.onUntrigger = onUntriggerFn;
      }
      window.tippy(el, config); 
    }
    const noterefs = window.document.querySelectorAll('a[role="doc-noteref"]');
    for (var i=0; i<noterefs.length; i++) {
      const ref = noterefs[i];
      tippyHover(ref, function() {
        // use id or data attribute instead here
        let href = ref.getAttribute('data-footnote-href') || ref.getAttribute('href');
        try { href = new URL(href).hash; } catch {}
        const id = href.replace(/^#\/?/, "");
        const note = window.document.getElementById(id);
        if (note) {
          return note.innerHTML;
        } else {
          return "";
        }
      });
    }
    const xrefs = window.document.querySelectorAll('a.quarto-xref');
    const processXRef = (id, note) => {
      // Strip column container classes
      const stripColumnClz = (el) => {
        el.classList.remove("page-full", "page-columns");
        if (el.children) {
          for (const child of el.children) {
            stripColumnClz(child);
          }
        }
      }
      stripColumnClz(note)
      if (id === null || id.startsWith('sec-')) {
        // Special case sections, only their first couple elements
        const container = document.createElement("div");
        if (note.children && note.children.length > 2) {
          container.appendChild(note.children[0].cloneNode(true));
          for (let i = 1; i < note.children.length; i++) {
            const child = note.children[i];
            if (child.tagName === "P" && child.innerText === "") {
              continue;
            } else {
              container.appendChild(child.cloneNode(true));
              break;
            }
          }
          if (window.Quarto?.typesetMath) {
            window.Quarto.typesetMath(container);
          }
          return container.innerHTML
        } else {
          if (window.Quarto?.typesetMath) {
            window.Quarto.typesetMath(note);
          }
          return note.innerHTML;
        }
      } else {
        // Remove any anchor links if they are present
        const anchorLink = note.querySelector('a.anchorjs-link');
        if (anchorLink) {
          anchorLink.remove();
        }
        if (window.Quarto?.typesetMath) {
          window.Quarto.typesetMath(note);
        }
        if (note.classList.contains("callout")) {
          return note.outerHTML;
        } else {
          return note.innerHTML;
        }
      }
    }
    for (var i=0; i<xrefs.length; i++) {
      const xref = xrefs[i];
      tippyHover(xref, undefined, function(instance) {
        instance.disable();
        let url = xref.getAttribute('href');
        let hash = undefined; 
        if (url.startsWith('#')) {
          hash = url;
        } else {
          try { hash = new URL(url).hash; } catch {}
        }
        if (hash) {
          const id = hash.replace(/^#\/?/, "");
          const note = window.document.getElementById(id);
          if (note !== null) {
            try {
              const html = processXRef(id, note.cloneNode(true));
              instance.setContent(html);
            } finally {
              instance.enable();
              instance.show();
            }
          } else {
            // See if we can fetch this
            fetch(url.split('#')[0])
            .then(res => res.text())
            .then(html => {
              const parser = new DOMParser();
              const htmlDoc = parser.parseFromString(html, "text/html");
              const note = htmlDoc.getElementById(id);
              if (note !== null) {
                const html = processXRef(id, note);
                instance.setContent(html);
              } 
            }).finally(() => {
              instance.enable();
              instance.show();
            });
          }
        } else {
          // See if we can fetch a full url (with no hash to target)
          // This is a special case and we should probably do some content thinning / targeting
          fetch(url)
          .then(res => res.text())
          .then(html => {
            const parser = new DOMParser();
            const htmlDoc = parser.parseFromString(html, "text/html");
            const note = htmlDoc.querySelector('main.content');
            if (note !== null) {
              // This should only happen for chapter cross references
              // (since there is no id in the URL)
              // remove the first header
              if (note.children.length > 0 && note.children[0].tagName === "HEADER") {
                note.children[0].remove();
              }
              const html = processXRef(null, note);
              instance.setContent(html);
            } 
          }).finally(() => {
            instance.enable();
            instance.show();
          });
        }
      }, function(instance) {
      });
    }
        let selectedAnnoteEl;
        const selectorForAnnotation = ( cell, annotation) => {
          let cellAttr = 'data-code-cell="' + cell + '"';
          let lineAttr = 'data-code-annotation="' +  annotation + '"';
          const selector = 'span[' + cellAttr + '][' + lineAttr + ']';
          return selector;
        }
        const selectCodeLines = (annoteEl) => {
          const doc = window.document;
          const targetCell = annoteEl.getAttribute("data-target-cell");
          const targetAnnotation = annoteEl.getAttribute("data-target-annotation");
          const annoteSpan = window.document.querySelector(selectorForAnnotation(targetCell, targetAnnotation));
          const lines = annoteSpan.getAttribute("data-code-lines").split(",");
          const lineIds = lines.map((line) => {
            return targetCell + "-" + line;
          })
          let top = null;
          let height = null;
          let parent = null;
          if (lineIds.length > 0) {
              //compute the position of the single el (top and bottom and make a div)
              const el = window.document.getElementById(lineIds[0]);
              top = el.offsetTop;
              height = el.offsetHeight;
              parent = el.parentElement.parentElement;
            if (lineIds.length > 1) {
              const lastEl = window.document.getElementById(lineIds[lineIds.length - 1]);
              const bottom = lastEl.offsetTop + lastEl.offsetHeight;
              height = bottom - top;
            }
            if (top !== null && height !== null && parent !== null) {
              // cook up a div (if necessary) and position it 
              let div = window.document.getElementById("code-annotation-line-highlight");
              if (div === null) {
                div = window.document.createElement("div");
                div.setAttribute("id", "code-annotation-line-highlight");
                div.style.position = 'absolute';
                parent.appendChild(div);
              }
              div.style.top = top - 2 + "px";
              div.style.height = height + 4 + "px";
              div.style.left = 0;
              let gutterDiv = window.document.getElementById("code-annotation-line-highlight-gutter");
              if (gutterDiv === null) {
                gutterDiv = window.document.createElement("div");
                gutterDiv.setAttribute("id", "code-annotation-line-highlight-gutter");
                gutterDiv.style.position = 'absolute';
                const codeCell = window.document.getElementById(targetCell);
                const gutter = codeCell.querySelector('.code-annotation-gutter');
                gutter.appendChild(gutterDiv);
              }
              gutterDiv.style.top = top - 2 + "px";
              gutterDiv.style.height = height + 4 + "px";
            }
            selectedAnnoteEl = annoteEl;
          }
        };
        const unselectCodeLines = () => {
          const elementsIds = ["code-annotation-line-highlight", "code-annotation-line-highlight-gutter"];
          elementsIds.forEach((elId) => {
            const div = window.document.getElementById(elId);
            if (div) {
              div.remove();
            }
          });
          selectedAnnoteEl = undefined;
        };
          // Handle positioning of the toggle
      window.addEventListener(
        "resize",
        throttle(() => {
          elRect = undefined;
          if (selectedAnnoteEl) {
            selectCodeLines(selectedAnnoteEl);
          }
        }, 10)
      );
      function throttle(fn, ms) {
      let throttle = false;
      let timer;
        return (...args) => {
          if(!throttle) { // first call gets through
              fn.apply(this, args);
              throttle = true;
          } else { // all the others get throttled
              if(timer) clearTimeout(timer); // cancel #2
              timer = setTimeout(() => {
                fn.apply(this, args);
                timer = throttle = false;
              }, ms);
          }
        };
      }
        // Attach click handler to the DT
        const annoteDls = window.document.querySelectorAll('dt[data-target-cell]');
        for (const annoteDlNode of annoteDls) {
          annoteDlNode.addEventListener('click', (event) => {
            const clickedEl = event.target;
            if (clickedEl !== selectedAnnoteEl) {
              unselectCodeLines();
              const activeEl = window.document.querySelector('dt[data-target-cell].code-annotation-active');
              if (activeEl) {
                activeEl.classList.remove('code-annotation-active');
              }
              selectCodeLines(clickedEl);
              clickedEl.classList.add('code-annotation-active');
            } else {
              // Unselect the line
              unselectCodeLines();
              clickedEl.classList.remove('code-annotation-active');
            }
          });
        }
    const findCites = (el) => {
      const parentEl = el.parentElement;
      if (parentEl) {
        const cites = parentEl.dataset.cites;
        if (cites) {
          return {
            el,
            cites: cites.split(' ')
          };
        } else {
          return findCites(el.parentElement)
        }
      } else {
        return undefined;
      }
    };
    var bibliorefs = window.document.querySelectorAll('a[role="doc-biblioref"]');
    for (var i=0; i<bibliorefs.length; i++) {
      const ref = bibliorefs[i];
      const citeInfo = findCites(ref);
      if (citeInfo) {
        tippyHover(citeInfo.el, function() {
          var popup = window.document.createElement('div');
          citeInfo.cites.forEach(function(cite) {
            var citeDiv = window.document.createElement('div');
            citeDiv.classList.add('hanging-indent');
            citeDiv.classList.add('csl-entry');
            var biblioDiv = window.document.getElementById('ref-' + cite);
            if (biblioDiv) {
              citeDiv.innerHTML = biblioDiv.innerHTML;
            }
            popup.appendChild(citeDiv);
          });
          return popup.innerHTML;
        });
      }
    }
  });
  </script>
</div> <!-- /content -->




</body></html><!--/html_preserve-->

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
[28] permute_0.9-8       htmltools_0.5.8.1   yaml_2.3.10        
[31] pillar_1.11.1       crayon_1.5.3        tidyr_1.3.1        
[34] MASS_7.3-65         vegan_2.7-2         iterators_1.0.14   
[37] foreach_1.5.2       nlme_3.1-168        tidyselect_1.2.1   
[40] digest_0.6.37       Rtsne_0.17          stringi_1.8.7      
[43] purrr_1.2.0         dplyr_1.1.4         reshape2_1.4.4     
[46] splines_4.5.2       fastmap_1.2.0       ade4_1.7-23        
[49] grid_4.5.2          cli_3.6.5           magrittr_2.0.4     
[52] survival_3.8-3      ape_5.8-1           withr_3.0.2        
[55] scales_1.4.0        XVector_0.50.0      igraph_2.2.1       
[58] multtest_2.66.0     evaluate_1.0.5      IRanges_2.44.0     
[61] mgcv_1.9-1          rlang_1.1.6         Rcpp_1.1.0         
[64] glue_1.8.0          BiocManager_1.30.26 renv_1.1.5         
[67] BiocGenerics_0.56.0 jsonlite_2.0.0      R6_2.6.1           
[70] Rhdf5lib_1.32.0     plyr_1.8.9         
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
