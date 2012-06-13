Pihur.RankMerging <-
function(exprSet,MergingDistance=c("Spearman", "Kendall"),...){
        PRLs=exprs(exprSet)
        MergingDistance<- match.arg(MergingDistance, c("Spearman", "Kendall"))
        # simple data preprocessing
        for (i in 1:ncol(PRLs))
           PRLs[,i]=as.matrix(rank(PRLs[,i]))

	phenodata=as(as(phenoData(exprSet),"data.frame"),"matrix")
	if (ncol(PRLs)!=length(phenodata))
		stop("the column of PRLs must be equal to the length of the phenodata")
	FPRL=matrix(0,nrow(PRLs))
	exp_names=phenodata
	if (ncol(phenodata)>1)
		phenodata=unique(phenodata,MARGIN=2)
	else
		phenodata=unique(phenodata,MARGIN=1)
	phenodata=sort(phenodata)
	phenodata_num=length(phenodata);
	gene_num=nrow(PRLs);
	exp_num=ncol(PRLs)
	tmp_indx=matrix()
	for (n1 in 1:phenodata_num){
		tmp_indx=matrix()
		diseasesI=phenodata[n1]
		k=1
		for (n2 in 1:length(exp_names)){
			if (diseasesI==exp_names[n2]){
				tmp_indx[k]=n2
				k=k+1
			}
		}
		R=PRLs[,tmp_indx];
		R=as.matrix(R)
                if (ncol(R)==1)
                        R=R
                else{
                       R=t(R)
                       if (gene_num>10){
                            if (gene_num<100)
                                R=RankAggreg(R,gene_num,weights=NULL,distance=MergingDistance,verbose=FALSE, ...)
                            else  stop("the size of rank aggregation is too large, please select the Iorio's method")
                       }
                       else  R=BruteAggreg(R,gene_num,weights=NULL,distance=MergingDistance, ...)
                       R=as.matrix(R$top.list) 
                 }
               FPRL=cbind(FPRL,R)
	       FPRL=as.matrix(FPRL)
       }
	FPRL=FPRL[,-1];
	phenodata=as.data.frame(phenodata)
        states=NULL
        for (i in 1:nrow(phenodata))
           states=cbind(states,paste(phenodata[i,1]))
         rownames(phenodata)=states
         colnames(FPRL)=rownames(phenodata)
	 phenodata=new("AnnotatedDataFrame",data=phenodata)
         return(new("ExpressionSet",exprs=FPRL,phenoData=phenodata))
}
