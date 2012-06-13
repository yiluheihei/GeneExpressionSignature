SignatureDistance <-
function(exprSet,SignatureLength,MergingMethod=c("Iorio","Pihur"),MergingDistance=c("Spearman", "Kendall"),
         ScoringMethod=c("GSEA","PGSEA"),ScoringDistance=c("avg","max"),...){
         MergingMethod=match.arg(MergingMethod,c("Iorio","Pihur"))
         MergingDistance=match.arg(MergingDistance, c("Spearman", "Kendall"))
         ScoringMethod=match.arg(ScoringMethod,c("GSEA","PGSEA"))
         ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
         if (MergingMethod=="Iorio")
             MergedRank=Iorio.RankMerging(exprSet,MergingDistance)
         else
             MergedRank=Pihur.RankMerging(exprSet,MergingDistance,...)
         if (ScoringMethod=="GSEA")
             return(ScoreGSEA(MergedRank,SignatureLength,ScoringDistance))
         else
             return(ScorePGSEA(MergedRank,SignatureLength,ScoringDistance))
}
