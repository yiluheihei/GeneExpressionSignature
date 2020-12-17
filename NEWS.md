# GeneExpressionSignature 1.37.0

* Add a `NEWS.md` file to track changes to the package.
* Add `inst/CITATION` file to customise the citation.
* Add `README.md`, `CODE_OF_CONDUCT.md`, and create a biocthis-style GitHub 
Actions workflow.
* Add a new function `PGSEA` from the [PGSEA](http://www.bioconductor.org/packages/PGSEA) to remove the dependency on 
[PGSEA](http://www.bioconductor.org/packages/PGSEA), since [PGSEA](http://www.bioconductor.org/packages/PGSEA) was deprecated in 
Bioconductor version 3.12 and removed from 3.13.
* Format code using [styler](https://cran.r-project.org/package=styler) 
and [biocthis](http://www.bioconductor.org/packages/biocthis), redocument package
using [Roxygen2](https://github.com/r-lib/roxygen2).
* Rewrite vignette using [knitr](), 
[rmarkdown](https://cran.r-project.org/package=styler), 
[BiocStyle](http://www.bioconductor.org/packages/BiocStyle).
* Fix BiocCheck errors and warnings.
