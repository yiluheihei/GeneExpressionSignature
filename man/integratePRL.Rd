\name{integratePRL}
\alias{integratePRL}
\title{Updating an existing dataset with new sample.
}
\description{Compute new enrichment score and distances among new PRL and previous PRLs in the existing dataset.
}
\usage{
integratePRL(ES, PRL, newPRL, qlen)
}
\arguments{
  \item{ES}{an ExpressionSet, array data is the existing enrichment score matrix
}
  \item{PRL}{the existing PRLs correspond to the ES
}
  \item{newPRL}{the PRL which you want to integrate to the existing PRLs 
}
  \item{qlen}{the length of the gene signature
}
}
\details{This function can integrate the new PRL into the previous PRLs to get the new enrichment score and distances matrix. 
}
\value{
\item{newPRLs }{an ExpressionSet, assay data is the PRL which new PRL have been integrated}
\item{newES }{an ExpressionSet, assay data is the integrated new ES matrix}
\item{newdistance }{an ExpressionSet, assay data is the integrated new distance matrix}
}
\seealso{ \code{\link{quickenrichmentscore}}
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

## compute distances and ES from aggregated matrix
d=distances(PRL,250)
ES=d[[1]]
distance=d[[2]]

## integrate new PRL to get newES and newdistances 
newPRL<- PRL[,2]
d <- integratePRL(ES,PRL,newPRL,250)
newES <- d[[2]]
newdistance <- d[[3]]

}


