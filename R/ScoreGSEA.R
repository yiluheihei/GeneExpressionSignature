ScoreGSEA <-
function(MergingSet,SignatureLength,DisMeasurement=c("avg", "max")){
        DisMeasurement=match.arg(DisMeasurement,c("avg","max"))
	ES=matrix(0) 
	FPRL=MergingSet
	phenodata=phenoData(MergingSet)
	if (is.null(ncol(exprs(FPRL)))==TRUE){
		stop("the class of argument FPRL is incorrect")
	}
	if (ncol(exprs(FPRL))<=1){
		stop("the number of column of argument FPRL must be greater than 1")
	}
	KBR=FPRL[,1]
	FPRLcol=ncol(exprs(FPRL))
	for (i in 2:FPRLcol){
		PRL=FPRL[,i]
		#PRL=as.matrix(PRL);
		d=integratePRL(ES,KBR,PRL,SignatureLength,DisMeasurement)
		ES=d[[2]]
		KBR=d[[1]]
	}
	return((d[[3]]))
        
}
