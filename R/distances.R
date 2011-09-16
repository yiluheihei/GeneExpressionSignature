distances <-
function(aggregateSet,qlen){
	ES=matrix(0) 
	FPRL=aggregateSet
	phenodata=phenoData(aggregateSet)
	ESpd=as(phenodata,"data.frame")
	ESpd=as.data.frame(ESpd[1,])
	ESpd=new("AnnotatedDataFrame",data=ESpd)
	ES=new("ExpressionSet",exprs=ES,phenoData=ESpd)
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
		d=integratePRL(ES,KBR,PRL,qlen)
		ES=d[[2]]
		KBR=d[[1]]
	}
	ES=d[[2]]
	distance=d[[3]]
	list(ES,distance)
}

