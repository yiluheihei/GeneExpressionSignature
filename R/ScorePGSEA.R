ScorePGSEA <-
function(MergingSet,SignatureLength,ScoringDistance=c("avg", "max"),p.value=F){
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
     Mvalue=max(abs(pgscores))
     pgscores=pgscores/max(abs(pgscores))  
     if (pgscores[1,1]>0){
       pgscores=pgscores
     }
     else{
       pgscores=-pgscores
     } 
     rownames(pgscores)=colnames(pgscores)
     p.results=matrix(0,n,n)
     if (ScoringDistance=="avg"){
        distances=1-(pgscores+t(pgscores))/2
        for (i in 1:n){
           for (j in 1:n){
               if (distances[i,j]<1)
                  p.results[i,j]=2 * pnorm(distances[i,j],mean=1,sd=1/(2*Mvalue))
               else
                  p.results[i,j]=2 * (1-pnorm(distances[i,j],mean=1,sd=1/(2*Mvalue)))
            }
        }
     }
     else{
       distances=pmin(1-pgscores,t(1-pgscores))/2
       p.results=matrix(0,nrow(distances),ncol(distances))
       for (i in 1:n){
           for (j in 1:n){
               if (distances[i,j]<1/2)
                  p.resutls[i,j]=2 * pnorm(distances[i,j],mean=1/2,sd=1/(8^(1/2)*Mvalue))
               else
                  p.resutls[i,j]=2 * (1-pnorm(distances[i,j],mean=1/2,sd=1/(8^(1/2)*Mvalue)))
            }
        }
     }
     if (p.value==T)
         return(list(distances = distances, p.results = p.results))
     return(distances)
}
