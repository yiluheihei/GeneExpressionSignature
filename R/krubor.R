krubor <-
function(distance=c("Spearman", "Kendall"),...){

        # find the two closest ranks and merge them into a rank
        findclosestrank <-function(SMDM)
        {
	     if (length(SMDM)!=0){
		   m=apply(SMDM,2,min);
		   i=apply(SMDM,2,which.min);
		   j=which.min(m);
	 	   m=min(m);
		   i=i[j];
 	      }
	     else{
		stop("zero-length SMDM is illegal")
	     }
	     list(m,i,j)
        }

       # Borda Mering method  
       BMRankMerging <-function(rankings)
       {
	   if (ncol(rankings)>1){
		majorities=rowSums(rankings);
		tmp=sort(majorities);
		sidxs=order(majorities);
		tmp=sort(sidxs);
		outrank=order(sidxs)
	    }
	   else
		outrank=rankings
	   return(as.matrix(outrank))
        }
      
      distance <- match.arg(distance, c("Spearman", "Kendall")) 
	R=data.frame(...)
	R=as.matrix(R)
	R=unique(R,MARGIN=2)
	R=as.matrix(R)
	nrank = ncol(R);
	while(nrank!=0) {
		SMDM=FootruleMatrix(R,distance,1)
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
