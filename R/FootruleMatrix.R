FootruleMatrix <-
function(Rankings,n){
	nrank=ncol(Rankings);
	if (length(Rankings)!=0){
		SMDM=matrix(0,nrow=nrank,ncol=nrank);
		if (nrank>1){
			for (i in 1:(nrank-1)) {
				for (j in (i+1):nrank){
					#optimize
					t=SMfootrule(Rankings[,i],Rankings[,j])
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

