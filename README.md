# GeneExpressionSignature: an R package for discovering functional connections using gene expression signatures

This package gives the implementations of the gene expression signature and 
its distance to each. Gene expression signature is represented as a list of 
genes whose expression is correlated with a biological state of interest. And 
its distance is defined using a nonparametric, rank-based pattern-matching 
strategy based on the Kolmogorov-Smirnov statistic. Gene expression signature 
and its distance can be used to detect similarities among the signatures of
drugs, diseases, and biological states of interest。

<!-- badges: start -->
<!-- badges: end -->


## Installation

You can install the released version of GeneExpressionSignature from [Bioconductor](http://www.bioconductor.org/packages/GeneExpressionSignature)
with:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(GeneExpressionSignature)
```

Or install the development version from github

```r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("yiluheihei/GeneExpressionSignature")
```

## Citation

**Authors**: Yang Cao <https://caoyang.tech>

If you use this package in your research, please cite:

Li F, Cao Y, Han L, Cui X, Xie D, Wang S, Bo X. [GeneExpressionSignature: an R package for discovering functional connections using gene expression 
signatures](https://doi.org/10.1089/omi.2012.0087). *Omics: a journal of integrative biology*, 2013 17(2):116-8.

## Contribution

Your contributions are always welcome!

Please note that the GeneExpressionSignature project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

