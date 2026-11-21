calc_in_frame<- function(normFactor, sample_names, file_list, output_location, s_rep_names){
  require(psych)
  require(reshape2)
  require(tidyverse)
  
  size_NF<- length(normFactor)
  print(file_list)
  normFactor_df <- data.frame(Replicate=s_rep_names, 
                              value_s=normFactor[grep("S", substring(names(normFactor),nchar(names(normFactor))))])
  if(length(grep("N", substring(names(normFactor),nchar(names(normFactor)))))>0){
    normFactor_df$value_n<- normFactor[grep("N", substring(names(normFactor),nchar(names(normFactor))))]
  }
  # normFactor_df <- data.frame(Replicate=s_rep_names, 
  #                             value_s=normFactor[grep("S", substring(names(normFactor),nchar(names(normFactor))))],
  #                             value_n=normFactor[grep("N", substring(names(normFactor),nchar(names(normFactor))))])
  # 
  # normFactor_df <- data.frame(Replicate=names(normFactor)[setdiff(c(1:size_NF),grep("N", names(normFactor)))], 
  #                             value_s=normFactor[setdiff(c(1:size_NF),grep("N", names(normFactor)))],
  #                             value_n=normFactor[grep("N", names(normFactor))])
  
  # filter_data <- function(count_table, norm_table){
  #   data <- count_table[count_table$num_fusion_reads_selected_sample>0 & count_table$num_fusion_reads_background_sample>0, 
  #                       c("Transcript_id", "Replicate", "num_fusion_reads_in_frame_selected_sample", 
  #                         "num_fusion_reads_selected_sample", "num_fusion_reads_in_frame_background_sample", "num_fusion_reads_background_sample")]
  #   data<- merge(data, norm_table, by="Replicate", all.x=T)
  #   data$num_fusion_reads_in_frame_selected<- data$num_fusion_reads_in_frame_selected_sample/data$value_s
  #   data$num_fusion_reads_selected_sample<- data$num_fusion_reads_selected_sample/data$value_s
  #   if("value_n" %in% colnames(data)){
  #     data$num_fusion_reads_in_frame_background<- data$num_fusion_reads_in_frame_background_sample/data$value_n
  #     data$num_fusion_reads_background<- data$num_fusion_reads_background/data$value_n
  #     data <- data %>% group_by(Transcript_id) %>%
  #       summarise(num_junction_reads_in_frame_s=round(sum(num_fusion_reads_in_frame_selected)),
  #                 num_junction_reads_s=round(sum(num_fusion_reads_selected_sample)),
  #                 num_junction_reads_in_frame_ns=round(sum(num_fusion_reads_in_frame_background)),
  #                 num_junction_reads_ns=round(sum(num_fusion_reads_background)))
  #     
  #     data$gene <- apply(data, 1, FUN=function(x){strsplit(x[1], "[.]")[[1]][1]})
  #     data <- data[data$num_junction_reads_s>=data$num_junction_reads_ns & data$num_junction_reads_ns>0,]
  #   }else{
  #     data <- data %>% group_by(Transcript_id) %>%
  #       summarise(num_junction_reads_in_frame_s=round(sum(num_fusion_reads_in_frame_selected)),
  #                 num_junction_reads_s=round(sum(num_fusion_reads_selected_sample)))
  #     
  #     data$gene <- apply(data, 1, FUN=function(x){strsplit(x[1], "[.]")[[1]][1]})
  #     data <- data[data$num_junction_reads_s>0,]
  #   }
  #   return(data)
  # }
  # 
  filter_data <- function(count_table, norm_table){
    if("Gene" %in% colnames(count_table)){
      data <- count_table[count_table$num_fusion_reads_selected_sample>0 & count_table$num_fusion_reads_background_sample>0, 
                          c("Gene","Transcript_id", "Replicate", "num_fusion_reads_in_frame_selected_sample", 
                            "num_fusion_reads_selected_sample", "num_fusion_reads_in_frame_background_sample", "num_fusion_reads_background_sample")]
    }else{
      data <- count_table[count_table$num_fusion_reads_selected_sample>0 & count_table$num_fusion_reads_background_sample>0, 
                          c("Transcript_id", "Replicate", "num_fusion_reads_in_frame_selected_sample", 
                            "num_fusion_reads_selected_sample", "num_fusion_reads_in_frame_background_sample", "num_fusion_reads_background_sample")]
      
    }
    
    data<- merge(data, norm_table, by="Replicate", all.x=T)
    data$num_fusion_reads_in_frame_selected<- data$num_fusion_reads_in_frame_selected_sample/data$value_s
    data$num_fusion_reads_selected_sample<- data$num_fusion_reads_selected_sample/data$value_s
    if(!"Gene" %in% colnames(data)){
      data$Gene <- apply(data, 1, FUN=function(x){strsplit(x[1], "[.]")[[1]][1]})
      }
    if("value_n" %in% colnames(data)){
      data$num_fusion_reads_in_frame_background<- data$num_fusion_reads_in_frame_background_sample/data$value_n
      data$num_fusion_reads_background<- data$num_fusion_reads_background/data$value_n
      data <- data %>% group_by(Transcript_id) %>%
        summarise(num_junction_reads_in_frame_s=round(sum(num_fusion_reads_in_frame_selected)),
                  num_junction_reads_s=round(sum(num_fusion_reads_selected_sample)),
                  num_junction_reads_in_frame_ns=round(sum(num_fusion_reads_in_frame_background)),
                  num_junction_reads_ns=round(sum(num_fusion_reads_background)),
                  gene=unique(Gene))
      
      
      data <- data[data$num_junction_reads_s>=data$num_junction_reads_ns & data$num_junction_reads_ns>0,]
    }else{
      data <- data %>% group_by(Transcript_id) %>%
        summarise(num_junction_reads_in_frame_s=round(sum(num_fusion_reads_in_frame_selected)),
                  num_junction_reads_s=round(sum(num_fusion_reads_selected_sample)),
                  gene=unique(Gene))
      
      #data$gene <- apply(data, 1, FUN=function(x){strsplit(x[1], "[.]")[[1]][1]})
      data <- data[data$num_junction_reads_s>0,]
    }
    return(data)
  }
  
  calculate_score <- function(count_table, alpha=0.05){
    data <- count_table
    data$in_frame_prop_s <- data$num_junction_reads_in_frame_s/data$num_junction_reads_s
    if("num_junction_reads_ns" %in% colnames(data)){
      data$in_frame_prop_ns <- data$num_junction_reads_in_frame_ns/data$num_junction_reads_ns
      data$in_frame_prop_ns_ho<- data$in_frame_prop_ns
      data$in_frame_prop_ns_ho[data$num_junction_reads_ns==0]<- 1/3
      data$pi_hat<- (data$in_frame_prop_s*data$num_junction_reads_s+data$in_frame_prop_ns_ho*data$num_junction_reads_ns)/
        (data$num_junction_reads_s+data$num_junction_reads_ns)
      data$statistic<- (data$in_frame_prop_s-data$in_frame_prop_ns_ho)/
        sqrt(data$pi_hat*(1-data$pi_hat)*(1/data$num_junction_reads_s+1/data$num_junction_reads_ns))
      
      
    }else{
      #assume 1/3 in ns samples, no info since no ns comtrol
      data$in_frame_prop_ns <- 1/3
      data$in_frame_prop_ns_ho<- data$in_frame_prop_ns
      data$statistic<- (data$in_frame_prop_s-data$in_frame_prop_ns_ho)/
        sqrt(data$in_frame_prop_s*(1-data$in_frame_prop_s)/data$num_junction_reads_s)
    }
    
    data<- data[!is.na(data$statistic),]
    #with normal
    #If we have data from NS libraries
    
    #data$p_val <- pnorm(q=data$statistic,lower.tail = FALSE)
    data$rank <- rank(data$statistic)
    data$freq_score <- data$rank/max(data$rank)
    data[is.na(data$freq_score),"freq_score"] <- 0
    return(data)
  }
  
  list_counts <- paste0(sample_names,"_final_report.csv")
  table_results_sum <- list()
  table_results <- list()
  for (num in 1:length(list_counts)){
    match <- sample_names[num]
    if(!(match %in% names(table))){
      t <- calculate_score(filter_data(read.csv(file_list[num],header = TRUE), normFactor_df))
      t_sum <- t %>% group_by(gene) %>% filter(num_junction_reads_in_frame_s> max(num_junction_reads_in_frame_s)/2) %>% summarise(total_max_freq_score=max(freq_score),
                                                  transcript= paste(Transcript_id[which(freq_score == max(freq_score))], collapse = ","))
      
      # t_sum <- t %>% group_by(gene) %>% summarise(total_max_freq_score=max(freq_score),
      #                                             transcript= paste(Transcript_id[which(freq_score == max(freq_score))], collapse = ","))
      
      table_results[[match]]<- t
      table_results[[match]]$bait<- match
      table_results_sum[[match]]<- t_sum
      table_results_sum[[match]]$bait<- match
      print(match)}
  } 
  in_frame_all_tab<-bind_rows(table_results)
  in_frame_results_sum<- bind_rows(table_results_sum)

  #saveRDS(in_frame_results_sum, paste0(output_location, "in_frame_results_sum.RDS"))
  return(in_frame_results_sum)
}





















