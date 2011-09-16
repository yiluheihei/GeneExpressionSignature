integratePRL <-
function(ES,PRL,newPRL,qlen){
	origpd=phenoData(PRL)
	newpd=phenoData(newPRL)
	ES=exprs(ES)
	PRL=exprs(PRL)
	prpd=new("AnnotatedDataFrame",data=rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
	newPRL=exprs(newPRL)
	PRL=cbind(PRL,newPRL)
	PRL=as.matrix(PRL)
	nelement=ncol(PRL)
	if (length(newPRL)>0) {
		n=length(newPRL);
		k1=1;k2=1
		UP=matrix()
		DOWN=matrix()
		for (j in 1:n){
			if (newPRL[j]<=qlen){
				UP[k1]=j
				k1=k1+1
			}
			if(newPRL[j]>=n-qlen+1){
				DOWN[k2]=j
				k2=k2+1
			}
		}
		ESrow=matrix(0,1,nelement)
		EScol=matrix(0,nelement,1)
		for (i in 1:nelement){
			ESrow[i]=quickenrichmentscore(UP,DOWN,PRL[,i])
			brgPRL=PRL[,i]
			brgPRL=as.matrix(brgPRL)
			n=nrow(brgPRL);
			k1=1;k2=1
			up=matrix()
			down=matrix()
			for (j in 1:n){
				if (brgPRL[j]<=qlen){
					up[k1]=j
					k1=k1+1
				}
				if(brgPRL[j]>=n-qlen+1){
					down[k2]=j
					k2=k2+1
				}
			}
			EScol[i]=quickenrichmentscore(up,down,newPRL)
		}
		ES=cbind(ES,EScol[c(1:(length(EScol)-1))])
		ES=rbind(ES,ESrow)
	}  
	ES=as.matrix(ES)
	distances=ES
	distances=as.matrix(distances)
	distances = (distances+t(distances))/2;
	distances = 1-distances;
	colnames(ES)=rownames(rbind(as(origpd,"data.frame"),as(newpd,"data.frame")))
	colnames(PRL)=colnames(ES)
	colnames(distances)=colnames(PRL)
	rownames(distances)=colnames(ES)
	rownames(ES)=colnames(ES)
	PRL=new("ExpressionSet",exprs=PRL,phenoData=prpd)
	ES=new("ExpressionSet",exprs=ES,phenoData=prpd)
	distances=new("ExpressionSet",exprs=distances,phenoData=prpd)
	list(PRL,ES,distances)
}

