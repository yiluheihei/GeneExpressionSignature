SignatureDistance <-
function(exprSet,SignatureLength,MergingDistance=c("Spearman", "Kendall"),
         ScoringMethod=c("GSEA","PGSEA"),ScoringDistance=c("avg","max"),...){
         MergingDistance=match.arg(MergingDistance, c("Spearman", "Kendall"))
         ScoringMethod=match.arg(ScoringMethod,c("GSEA","PGSEA"))
         ScoringDistance=match.arg(ScoringDistance,c("avg","max"))
         MergedRank=RankMerging(exprSet,MergingDistance)
         if (ScoringMethod=="GSEA")
             return(ScoreGSEA(MergedRank,SignatureLength,ScoringDistance))
         else
             return(ScorePGSEA(MergedRank,SignatureLength,ScoringDistance))
}
