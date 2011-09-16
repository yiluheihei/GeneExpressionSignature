\name{aggregate}
\alias{aggregate}
\title{
Aggregate each group of ranked lists with the same state into a single list
}
\description{
Aggregate the assay data according to phenotypic data of the input ExpressionSet.  Each group of the ranked lists with the same phenotypic data is aggregated into a single list, return it as an ExpressionSet object.
}
\usage{
aggregate(exprSet)
}
\arguments{
	\item{exprSet}{an ExpressionSet object, each column of assay data represents a ranked list obtained by preprocessing the corresponding gene expression profile, and phenotypic data represents the short description (characteristics of gene expression profile, such as the drug type, the disease state) about the assay data. 
	}
}
\details{
The krubor function is used in the aggregating procedure. And the following methods are used in the implementation: a measure of the distance between two ranked lists (Spearman's Footrule), a method to merge two or more ranked lists the (Borda Merging Method), and a algorithm to obtain a single ranked list from a set of them in a hierarchical way (the Kruskal Algorithm).
}

\seealso{\code{\link{krubor},aggregate all ranked lists into one list }}

\examples{
library(Biobase)
## load sample ranked list
PRLs=as.matrix(read.table(system.file("extdata/example_PRLs.txt",package="GeneExpressionSignature")))
## load sample phenotypic data
states=read.table(system.file("extdata/example_states.txt",package="GeneExpressionSignature"))
## create an new ExpressionSet object
rownames(states)=colnames(PRLs)
phenodata=new("AnnotatedDataFrame",data = states)
exprSet=new("ExpressionSet",exprs=PRLs,phenoData=phenodata)
## aggregate each group of the ranked lists in the exprSet with the same phenotypic data into a single PRL
aggregateSet=aggregate(exprSet)
}

