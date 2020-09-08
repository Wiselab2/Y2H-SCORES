calculate_enrichment_score<- function(dds, sample_names, p_val_thr, output_location){
  require(tidyverse)
  contrast_tables<- list()
  enrichment_score<- data.frame()
  for(sample in sample_names){
    contrast<-data.frame(results(dds, contrast=c("group", paste0(sample,"S"), paste0(sample,"N"))),
                         cooksCutoff=FALSE)
    contrast_tables[[sample]]<- contrast[!is.na(contrast$pvalue),]
    contrast_tables[[sample]]$gene<- rownames(contrast_tables[[sample]])
    contrast_tables[[sample]]$bait<- sample
  }
  
  #now we get the information to build the score
  bind_contrast_tables<- bind_rows(contrast_tables)
  enrichment_score<- bind_contrast_tables[bind_contrast_tables$pvalue<p_val_thr & bind_contrast_tables$log2FoldChange>0,]
  enrichment_score$rank <- rank(-enrichment_score$stat)
  #here for new thresholds
  enrichment_score$bait_enrich_score <- (max(enrichment_score$rank)-enrichment_score$rank)/max(enrichment_score$rank)
  #saveRDS(enrichment_score, paste0(output_location,"enrichment_score.RDS"))
  #saveRDS(bind_contrast_tables, paste0(output_location,"contrast_tables_enrichment_score.RDS"))
  get_rel_fc_score <- function(score_table, x, y, n = 100) {
    library(MASS)
    dens <- MASS::kde2d(x = score_table[x][[1]], y = score_table[y][[1]], n = n, h= rep(0.6,2))
    table <- score_table
    table[paste0("rel_",y,"_score")] <-0
    table$index <- findInterval(score_table[x][[1]], dens$x)
    for (j in 1:n) {
      table[table$index==j,paste0("rel_",y,"_score")]<- rank(table[table$index==j,y])/max(rank(table[table$index==j,y]))
      table[table$index==j,paste0("total_",x)]<- table[table$index==j,x] +
        table[table$index==j,paste0("rel_",y,"_score")]*(max(table[table$index==j,x])-min(table[table$index==j,x])+0.0001)/sum(table$index==j)}
    return(table)
  }
  
  enrichment_score_all_rel <- get_rel_fc_score(enrichment_score, "bait_enrich_score", "log2FoldChange")
  #saveRDS(enrichment_score_all_rel, paste0(output_location,"enrichment_score_all_rel.RDS"))
  return(enrichment_score_all_rel)
}

