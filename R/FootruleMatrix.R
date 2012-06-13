FootruleMatrix <-
function(Rankings,MergingDistance=c("Spearman", "Kendall"),n){

        # Kendall distance,this source code is derived from R package RankAggreg
        Kendallfootrule <-function(x, y)
        {
           K=0
           n <- length(x)
           for (i in 1:(n-1))
               for (j in i:n)
                   if((x[i] > x[j] & y[i] < y[j]) | (x[i] < x[j] & y[i] > y[j]))
                      K=K+1
           K
        } 
       
       # Spearman distance   
         SMfootrule <-function(R1,R2)
         {
	   SMD=sum(abs(R1-R2));
	   SMD
         }

        MergingDistance <- match.arg(MergingDistance, c("Spearman", "Kendall"))
	nrank=ncol(Rankings);
	if (length(Rankings)!=0){
		SMDM=matrix(0,nrow=nrank,ncol=nrank);
		if (nrank>1){
			for (i in 1:(nrank-1)) {
				for (j in (i+1):nrank){
					#optimize
                                        if(MergingDistance=="Spearman")
					   t=SMfootrule(Rankings[,i],Rankings[,j])
                                        else
                                           t=Kendallfootrule(Rankings[,i],Rankings[,j])
					if (t>0)
						SMDM[i,j]=t
					else
						SMDM[i,j]=0
				}
			}
			SMDM=SMDM+t(SMDM);
			if (n!=0){
				SMDMrol=nrow(SMDM);
				SMDMcol=ncol(SMDM);
				if (max(SMDMrol,SMDMcol)>1){
					if (max(SMDM[])>0)
						SMDM=SMDM/max(SMDM[])
				}
			}
		}
		else if (nrank==1){
			SMDM=1
		}
		else
			SMDM=0
	}
	return(SMDM)
}
