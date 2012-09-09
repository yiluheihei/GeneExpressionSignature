\name{getRLs}
\alias{getRLs}
\title{
Convert gene expression profiles to a ranked list (data preprocessing)
}
\description{
Sorting the microarray probe-set identifiers according to the differential expression values with respect to the untreated hybridization to obtaine a ranked list. Gene-expression profiles in are represented in a nonparametric fashion. 
}
\usage{
getRLs(control, treatment)
}
\arguments{
  \item{control}{
a matrix, including the vehicle control gene expression profiles correspongding to the treatment gene expression profiles. 
}
  \item{treatment}{
a matrix, is composed of gene expression profiles. 
}
}
\details{
The genes on the array are rank-ordered according to their differential expression relative to the control. First, control and treatment values less than a primary threshold value (quartile) were set to that threshold value. Finally, probe sets were ranked in descending order of d, where d is the ratio of the corresponding treatment-to-control values. For probe sets where d=1, a lower threshold was applied to the original difference values and a new treatment to control ratio (d') calculated. These probe sets were then sub-sorted in descending order of d.
}
\value{
A matrix is composed of ranked lists, a ranked list represents the correspongding gene expression profiles.
}
\examples{

if (require(GEOquery)){
  #treatment gene-expression profiles
  GSM118720 <- getGEO(filename=system.file("extdata/GSM118720.soft",package=
  "GeneExpressionSignature"))
  #control gene-expression profiles
  GSM118721 <- getGEO(filename=system.file("extdata/GSM118721.soft",package=
  "GeneExpressionSignature"))
  #data ranking according to the different expression values 
  control <- as.matrix(as.numeric(Table(GSM118721)[,2]))
  treatment <- as.matrix(as.numeric(Table(GSM118720)[,2]))
  ranked_list <-getRLs(control,treatment)
}

}
