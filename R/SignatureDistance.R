SignatureDistance <-
function(exprSet,MergingMethod=c("Iorio","Pihur"),distance=c("Spearman", "Kendall"),
         ScoringMethod=c("GSEA","PGSEA"),SignatureLength,DisMeasurement=c("avg","max"),...){
         MergingMethod=match.arg(MergingMethod,c("Iorio","Pihur"))
         distance=match.arg(distance, c("Spearman", "Kendall"))
         ScoringMethod=match.arg(ScoringMethod,c("GSEA","PGSEA"))
         DisMeasurement=match.arg(DisMeasurement,c("avg","max"))
         if (MergingMethod=="Iorio")
             MergedRank=Iorio.RankMerging(exprSet,distance)
         else
             MergedRank=Pihur.RankMerging(exprSet,distance,...)
         if (ScoringMethod=="GSEA")
             return(ScoreGSEA(MergedRank,SignatureLength, DisMeasurement))
         else
             return(ScorePGSEA(MergedRank,SignatureLength, DisMeasurement))
}
