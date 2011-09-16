aggregate <-
function(exprSet){
	PRLs=exprs(exprSet)
	phenodata=as(as(phenoData(exprSet),"data.frame"),"matrix")
	if (ncol(PRLs)!=length(phenodata))
		stop("the column of PRLS must be equal to the length of the phenodata")
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
		R=krubor(R)
		FPRL=cbind(FPRL,R)
		FPRL=as.matrix(FPRL)
	}
	FPRL=FPRL[,-1];
	phenodata=as.data.frame(phenodata)
	colnames(FPRL)=c(1:ncol(FPRL))
	rownames(phenodata)=c(1:nrow(phenodata))
	phenodata=new("AnnotatedDataFrame",data=phenodata)
	return(new("ExpressionSet",exprs=FPRL,phenoData=phenodata))
}

