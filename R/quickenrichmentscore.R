quickenrichmentscore <-
function(S,S1,List){
	Or_list=List;
	#List=sort(List);
	#List=as.matrix(List)
	Rank=order(Or_list)
	Rank=as.matrix(Rank)
	N=length(Rank)
	Nh=length(S)
	tmp=matrix(0,N)
	for (i in 1:Nh){
		tmp[Or_list[S[i]]]=1
	}
	hitCases=cumsum(tmp)
	missCases=cumsum(1-tmp)
	NR=length(S)
	Phit=hitCases/NR
	Pmiss=missCases/(N-Nh)
	abshm=abs(Phit-Pmiss)
	abshm=as.matrix(abshm)
	#m=apply(abshm,2,max);
	t=apply(abshm,2,which.max)
	ESUP = Phit[t]-Pmiss[t];
	RS = Phit-Pmiss;

	Or_list2=Or_list
	Rank=order(Or_list2)
	Rank=as.matrix(Rank) 
	N=length(Rank)
	Nh=length(S1)
	tmp=matrix(0,N)
	for (i in 1:Nh){
		tmp[Or_list[S1[i]]]=1
	}
	hitCases=cumsum(tmp)
	missCases=cumsum(1-tmp) 
	NR=length(S1)
	Phit=hitCases/NR
	Pmiss=missCases/(N-Nh)
	abshm=abs(Phit-Pmiss)
	abshm=as.matrix(abshm) 
	#m=apply(abshm,2,max);
	t=apply(abshm,2,which.max)
	ESDOWN = Phit[t]-Pmiss[t];
	RS = Phit-Pmiss;
	ES = (ESUP - ESDOWN)/2;
}
