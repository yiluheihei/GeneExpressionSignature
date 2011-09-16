BMRankMerging <-
function(rankings){
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

