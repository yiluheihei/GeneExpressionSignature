krubor <-
function(...){
	R=data.frame(...)
	R=as.matrix(R)
	R=unique(R,MARGIN=2)
	R=as.matrix(R)
	nrank = ncol(R);
	while(nrank!=0) {
		SMDM=FootruleMatrix(R,1)
		SMDM=as.matrix(SMDM)
		SMDM[lower.tri(SMDM)]=0
		if (nrank==1)
			nrank=0
		else{
			SMDM[SMDM==0]='inf'
			mij=findclosestrank(SMDM)
			i=mij[[2]]
			j=mij[[3]]
			R1=R[,i];
			R2=R[,j];
			R12=cbind(R1,R2);
			R12=as.matrix(R12)
			newrank= BMRankMerging(R12);
			newrank=as.matrix(newrank)
			nrank =nrank-1;
			R=R[,c(-i,-j)];
			R=cbind(R,newrank);
		}
	}
	R
}

