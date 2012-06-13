integratePRL <-
function(ES,PRL,newPRL,SignatureLength,ScoringDistance=c("avg", "max")){
        ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
	origpd=phenoData(PRL)
	newpd=phenoData(newPRL)
	#ES=exprs(ES)
	PRL=exprs(PRL)
	prpd=new("AnnotatedDataFrame",data=rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
	newPRL=exprs(newPRL)
	PRL=cbind(PRL,newPRL)
	PRL=as.matrix(PRL)
	nelement=ncol(PRL)
	if (length(newPRL)>0) {
		n=length(newPRL);
                UP=which(newPRL<=SignatureLength)
                DOWN=which(newPRL>=n-SignatureLength+1)
		ESrow=matrix(0,1,nelement)
		EScol=matrix(0,nelement,1)
		for (i in 1:nelement){
			ESrow[i]=quickenrichmentscore(UP,DOWN,PRL[,i])
			brgPRL=PRL[,i]
			brgPRL=as.matrix(brgPRL)
			n=nrow(brgPRL);
                        up=which(brgPRL<=SignatureLength)
                        down=which(brgPRL>=n-SignatureLength+1)
			EScol[i]=quickenrichmentscore(up,down,newPRL)
		}
		ES=cbind(ES,EScol[c(1:(length(EScol)-1))])
		ES=rbind(ES,ESrow)
	}  
	ES=as.matrix(ES)
	distances=ES
	distances=as.matrix(distances)
        if (ScoringDistance=="avg")
              distances = (distances+t(distances))/2
        else
              distances = pmax(distances+t(distances))
	distances = 1-distances
	colnames(ES)=rownames(rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
	colnames(PRL)=colnames(ES)
	colnames(distances)=colnames(PRL)
	rownames(distances)=colnames(ES)
	rownames(ES)=colnames(ES)
	PRL=new("ExpressionSet",exprs=PRL,phenoData=prpd)
	list(PRL,ES,distances)
}
