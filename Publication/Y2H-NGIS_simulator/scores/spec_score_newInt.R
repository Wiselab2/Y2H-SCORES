
calc_spec_score<- function(dds, sample_names, p_val_thr_s, fc_thr, output_location){
  require(DESeq2)
  require(tidyverse)
  combi <-combn(sample_names,2)
  spec_tables<- list()
  for(i in 1:dim(combi)[2]){
    contrast<-data.frame(results(dds, contrast=c("group", paste0(combi[1,i],"S"),
                                                 paste0(combi[2,i],"S")), cooksCutoff = FALSE))
    contrast$gene <- rownames(contrast)
    contrast$bait_num <- combi[1,i]
    contrast$bait_den <- combi[2,i]
    
    if(is.null(nrow(spec_tables[[combi[1, i]]]))){
      spec_tables[[combi[1,i]]]<- contrast[contrast$log2FoldChange>0,]
    }else{
      spec_tables[[combi[1,i]]]<- rbind(spec_tables[[combi[1,i]]],
                                        contrast[contrast$log2FoldChange>0,])
    }
    if(is.null(nrow(spec_tables[[combi[2, i]]]))){
      spec_tables[[combi[2,i]]]<- contrast[contrast$log2FoldChange<0,]
    }else{
      spec_tables[[combi[2,i]]]<- rbind(spec_tables[[combi[2,i]]],
                                        contrast[contrast$log2FoldChange<0,])
    }
    print(i)
  }
  for (sample in sample_names) {
    spec_tables[[sample]]$bait<-sample
  }
  
  spec_tables_all <- bind_rows(spec_tables)
  spec_tables_all<- spec_tables_all[!is.na(spec_tables_all$gene)&!duplicated(spec_tables_all),]
  saveRDS(spec_tables_all, paste0(output_location,"spec_tables_all.RDS"))
  spec_score<- spec_tables_all[spec_tables_all$pvalue<p_val_thr_s & abs(spec_tables_all$log2FoldChange)>fc_thr,]
  #now we get the information to build the score
  #spec_score$rank <- rank(spec_score$pvalue)
  spec_score$rank <- rank(-abs(spec_score$stat))
  
  #here for new thresholds
  spec_score$bait_spec_score <- (max(spec_score$rank)-spec_score$rank)/max(spec_score$rank)
  
  get_rel_fc_score2 <- function(score_table, x, y, n = 100) {
    library(MASS)
    dens <- MASS::kde2d(x = score_table[x][[1]], y = abs(score_table[y][[1]]), n = n, h= rep(0.6,2))
    table <- score_table
    table[paste0("rel_",y,"_score")] <-0
    table$index <- findInterval(score_table[x][[1]], dens$x)
    for (j in 1:n) {
      table[table$index==j,paste0("rel_",y,"_score")]<- rank(table[table$index==j,y])/
        max(rank(table[table$index==j,y]))
      table[table$index==j,paste0("rel_",y,"_score_interval")]<-table[table$index==j,paste0("rel_",y,"_score")]*
        (max(table[table$index==j,x])- min(table[table$index==j,x])+0.01)/sum(table$index==j)
      table[table$index==j,paste0("subtotal_",x)]<- table[table$index==j,x] +
        table[table$index==j,paste0("rel_",y,"_score")]*(max(table[table$index==j,x])-
                                                           min(table[table$index==j,x])+0.0001)/sum(table$index==j)}
    return(table)
  }
  
  spec_score_comp2<- data.frame(spec_score %>% group_by(gene,bait) %>% mutate(spec_score= sum(bait_spec_score)/length(sample_names), log2FoldChange_abs=abs(log2FoldChange)))
  spec_score_comp <- get_rel_fc_score2(spec_score_comp2, "spec_score","log2FoldChange_abs")
  spec_score_total<- spec_score_comp %>% group_by(gene,bait) %>% summarise(total_spec_score=mean(spec_score)+ sum(rel_log2FoldChange_abs_score_interval)/length(sample_names))
  saveRDS(spec_score_total, paste0(output_location,"spec_score_all_rel.RDS"))
  return(spec_score_total)
}
