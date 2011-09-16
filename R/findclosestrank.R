findclosestrank <-
function(SMDM){
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

