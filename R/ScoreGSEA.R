ScoreGSEA <-
function(MergingSet,SignatureLength,ScoringDistance=c("avg", "max")){
        ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
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
		d=integratePRL(ES,KBR,PRL,SignatureLength,ScoringDistance)
		ES=d[[2]]
		KBR=d[[1]]
	}
	return((d[[3]]))
        
}
