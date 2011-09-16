\name{distances}
\alias{distances}
\title{
Compute pairwise distances between samples
}
\description{
Compute pairwise distances between sample according to their (Prototype Ranked List) PRL, get a n-n distance matrix as the assay data of the result , n
is the length of PRL. 
}
\usage{
distances(aggregateSet, qlen)
}
\arguments{
  \item{aggregateSet}{an ExpressionSet object. The assay data represents the PRLs of the samples, each column represents one PRL. The number of column of this argument must be greater than
1, oherwise, this function is not meaningful.
}
  \item{qlen}{the length of "gene signature". In order to compute pairwise distances among samples, genes lists are ranked according to the gene expression ratio (fold change). And the "gene signature" includes the most up-regulated genes (near the top of the list) and the most down-regulated genes (near the bottom of the list). 
}
}
\details{Once the PRL obtained for each sample, the distances between samples are calculated base on gene signature, including the expression of genes that seemed to consistently vary in response to the across different experimental conditions (e.g., different cell lines and different dosages).
}
\value{
\item{ES }{an ExpressionSet, assay data is the enrichment score matrix}
\item{DS }{an ExpressionSet, assay data is the distance matrix, the maximum distance is more sensitive to weak similarities, providing a lower precision but a larger recall.}
}
\seealso{\code{\link{aggregate}}
}
\examples{
## create an instance ExpressionSet
library(Biobase)
PRLs=as.matrix(read.table(system.file("extdata/example_PRLs.txt",package="GeneExpressionSignature")))
states=read.table(system.file("extdata/example_states.txt",package="GeneExpressionSignature"))
rownames(states)=colnames(PRLs)
phenodata=new("AnnotatedDataFrame",data = states)
exampleSet=new("ExpressionSet",exprs=PRLs,phenoData=phenodata)

## aggregate the exampleSet 
PRL=aggregate(exampleSet)

## compute distances from aggregated matrix
d=distances(PRL,250)
enrichmentscore=d[[1]]
distance=d[[2]]
}


