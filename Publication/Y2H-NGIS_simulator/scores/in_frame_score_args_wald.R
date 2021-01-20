calc_in_frame<- function(normFactor, sample_names, file_list, output_location){
  require(psych)
  require(reshape2)
  require(tidyverse)
  
  size_NF<- length(normFactor)
  print(file_list)
  normFactor_df <- data.frame(Replicate=names(normFactor)[setdiff(c(1:size_NF),grep("N", names(normFactor)))], 
                              value_s=normFactor[setdiff(c(1:size_NF),grep("N", names(normFactor)))],
                              value_n=normFactor[grep("N", names(normFactor))])
  
  filter_data <- function(count_table, norm_table){
    data <- count_table[count_table$num_fusion_reads_selected_sample>0 & count_table$num_fusion_reads_background_sample>0, 
                        c("Transcript_id", "Replicate", "num_fusion_reads_in_frame_selected_sample", 
                          "num_fusion_reads_selected_sample", "num_fusion_reads_in_frame_background_sample", "num_fusion_reads_background_sample")]
    data<- merge(data, norm_table, by="Replicate", all.x=T)
    data$num_fusion_reads_in_frame_selected<- data$num_fusion_reads_in_frame_selected_sample/data$value_s
    data$num_fusion_reads_selected_sample<- data$num_fusion_reads_selected_sample/data$value_s
    data$num_fusion_reads_in_frame_background<- data$num_fusion_reads_in_frame_background_sample/data$value_n
    data$num_fusion_reads_background<- data$num_fusion_reads_background/data$value_n
    
    data <- data %>% group_by(Transcript_id) %>%
      summarise(num_junction_reads_in_frame_s=round(sum(num_fusion_reads_in_frame_selected)),
                num_junction_reads_s=round(sum(num_fusion_reads_selected_sample)),
                num_junction_reads_in_frame_ns=round(sum(num_fusion_reads_in_frame_background)),
                num_junction_reads_ns=round(sum(num_fusion_reads_background)))
    
    data$gene <- apply(data, 1, FUN=function(x){strsplit(x[1], "[.]")[[1]][1]})
    data <- data[data$num_junction_reads_s>=data$num_junction_reads_ns & data$num_junction_reads_ns>0,]
    return(data)
  }
  
  calculate_score <- function(count_table, alpha=0.05){
    data <- count_table
    data$in_frame_prop_s <- data$num_junction_reads_in_frame_s/data$num_junction_reads_s
    data$in_frame_prop_ns <- data$num_junction_reads_in_frame_ns/data$num_junction_reads_ns
    
    #with normal
    #If we have data from NS libraries
    data$in_frame_prop_ns_ho<- data$in_frame_prop_ns
    data$in_frame_prop_ns_ho[data$num_junction_reads_ns==0]<- 1/3
    data$statistic<- (data$in_frame_prop_s-data$in_frame_prop_ns_ho)/
      sqrt(data$in_frame_prop_s*(1-data$in_frame_prop_s)/(data$num_junction_reads_s)+
             data$in_frame_prop_ns_ho*(1-data$in_frame_prop_ns_ho)/(data$num_junction_reads_ns))
    data<- data[!is.na(data$statistic),]
    data$p_val <- pnorm(q=data$statistic,lower.tail = FALSE)
    data$rank <- rank(data$statistic)
    data$freq_score <- data$rank/max(data$rank)
    #data$freq_score <- pnorm(q=data$statistic) #1-data$p_val
    data$plog2FC <- (data$num_junction_reads_in_frame_s+1)/(data$num_junction_reads_in_frame_ns+1)
    
    data$statistic_np<- (data$num_junction_reads_in_frame_s-data$num_junction_reads_in_frame_ns)/
      sqrt(1+data$in_frame_prop_s*(1-data$in_frame_prop_s)*(data$num_junction_reads_s)+
             data$in_frame_prop_ns_ho*(1-data$in_frame_prop_ns_ho)*(data$num_junction_reads_ns))
    data<- data[!is.na(data$statistic_np),]
    data$rank_np <- rank(data$statistic_np)
    data$freq_score_np <- data$rank_np/max(data$rank_np)
    data[is.na(data$freq_score),"freq_score"] <- 0
    data[is.na(data$freq_score_np),"freq_score_np"] <- 0
    return(data)
  }
  
  get_rel_IF_score <- function(score_table, x, y, n = 100) {
    library(MASS)
    dens <- MASS::kde2d(x = score_table[x][[1]], y = score_table[y][[1]], n = n, h= rep(0.6,2))
    table <- score_table
    table[paste0("rel_",y,"_score")] <-0
    table$index <- findInterval(score_table[x][[1]], dens$x)
    for (j in unique(table$index)) {
      table[table$index==j,paste0("rel_",y,"_score")]<- rank(table[table$index==j,y][[1]])/
        max(rank(table[table$index==j,y][[1]]))
      table[table$index==j,paste0("total_",x)]<- table[table$index==j,x] +
        table[table$index==j,paste0("rel_",y,"_score")]*(max(table[table$index==j,x])- min(table[table$index==j,x])+0.01)}
    table[paste0("total_",x)]<-table[paste0("total_",x)]/max(table[paste0("total_",x)])
    return(table)
  }
  
  list_counts <- paste0(sample_names,"_final_report.csv")
  table_results_sum <- list()
  table_results <- list()
  for (num in 1:length(list_counts)){
    match <- sample_names[num]
    if(!(match %in% names(table))){
      t <- calculate_score(filter_data(read.csv(file_list[num],header = TRUE), normFactor_df))
      t_sum <- t %>% group_by(gene) %>% summarise(max_freq_score_np=max(freq_score_np[which(plog2FC == max(plog2FC))]), 
                                                  max_freq_score=max(freq_score[which(plog2FC == max(plog2FC))]))
      table_results[[match]]<- t
      table_results[[match]]$bait<- match
      table_results_sum[[match]]<- t_sum
      table_results_sum[[match]]$bait<- match
      print(match)}
  } 
  in_frame_all_tab<-bind_rows(table_results)
  in_frame_results_sum<- bind_rows(table_results_sum)
  in_frame_results_sum_np_freq <- get_rel_IF_score(in_frame_results_sum, "max_freq_score_np", "max_freq_score")
  
  saveRDS(in_frame_results_sum_np_freq, paste0(output_location, "in_frame_results_sum_rel_np_freq.RDS"))
  return(in_frame_results_sum_np_freq)
  }





















