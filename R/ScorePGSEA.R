ScorePGSEA <-
function(MergingSet,SignatureLength,ScoringDistance=c("avg", "max")){
    ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
    PRLs=exprs(MergingSet)
    n=ncol(PRLs)
    m=nrow(PRLs)
    pgscores=matrix(0,1,n)
    for (i in 1:n){
        prls=as.matrix(PRLs[,i])
        up=which(prls<=SignatureLength)
        down=which(prls>=m-SignatureLength+1)
        upsmc=new("smc",ids=up)
        downsmc=new("smc",ids=down)
        uppgscore=PGSEA(MergingSet,list(upsmc),p.value=NA)
        downpgscore=PGSEA(MergingSet,list(downsmc),p.value=NA)
        pgscore=(uppgscore-downpgscore)/2
        pgscores=rbind(pgscores,pgscore)
        
    }
     pgscores=as.matrix(pgscores)
     pgscores=pgscores[-1,]
     pgscores=pgscores/max(abs(pgscores))  
     if (pgscores[1,1]>0){
       pgscores=pgscores
     }
     else{
       pgscores=-pgscores
     } 
     rownames(pgscores)=colnames(pgscores)
     if (ScoringDistance=="avg")
       distances=1-(pgscores+t(pgscores))/2
     else
       distances=pmax(1-pgscores,t(1-pgscores))
    return(distances)   
}
