run_scores<-function(sample_names, num_replicates, input_files_location, location_folders, out_dir, p_val_thr, p_val_thr_s, fc_thr, norm_plots=F){
  require(psych)
  require(reshape2)
  require(DESeq2)
  salmon_counts_list<-paste(input_files_location,sample_names,"_salmon_counts.matrix",sep="")
  final_reports_list<- paste(input_files_location,sample_names,"_final_report.csv",sep="")
  dir.create(out_dir)
  print(p_val_thr_s)
  print(fc_thr)
  print(p_val_thr)
  print(out_dir)
  print(salmon_counts_list)
  print(final_reports_list)
  print(sample_names)
  print(num_replicates)
  print(str(num_replicates))
  print(location_folders)
  
  ############################################################################
  output_location<- out_dir
  # Read in the data from the files and create a data frame for the entire data
  raw_counts<-data.frame()
  for(sample_num in 1:length(sample_names)){
    raw_counts_per_sample<-read.table(salmon_counts_list[sample_num], row.names = 1)
    column_names<-c(paste0("R",1:num_replicates[sample_num], "S"), paste0("R",1:num_replicates[sample_num], "N"))
    colnames(raw_counts_per_sample)<-paste(sample_names[sample_num],column_names,sep="")
    if(nrow(raw_counts)==0)
    {
      raw_counts<-raw_counts_per_sample
    }else{
      raw_counts<-cbind(raw_counts,raw_counts_per_sample)
    }
  }
  rm(raw_counts_per_sample)
  
  # Remove all those genes which have a count of 0 across ALL samples
  raw_counts<-raw_counts[rowSums(raw_counts)>0,]
  # Convert counts to integers
  raw_counts<-round(raw_counts)
  totCounts <- colSums(raw_counts)
  normFactor <- totCounts/geometric.mean(totCounts)
  normalized_counts<- raw_counts/rep(normFactor, each = (nrow(raw_counts)))
  normalized_counts<-round(normalized_counts)
  saveRDS(normFactor, paste0(output_location,"normFactor_salmon.RDS"))
  saveRDS(raw_counts, paste0(output_location,"raw_counts_salmon.RDS"))
  #run in frame
  source(paste0(location_folders,"/in_frame_score_args_wald.R"))
  in_frame_all<- calc_in_frame(normFactor, sample_names, final_reports_list, output_location)
  # Run DESeq2 for enriched and specificity score
  condition<-c("S","N")
  cols<-data.frame(baits=rep(sample_names, 2*num_replicates), 
                   conditions=unlist(list(apply(data.frame(num_replicates), 1, FUN=function(x){rep(condition,each=x[1])}))), 
                   replication=unlist(list(apply(data.frame(num_replicates), 1, FUN=function(x){rep(1:x[1],2)}))))
  rownames(cols)<-colnames(raw_counts)
  cols$conditions<-as.factor(cols$conditions)
  cols$replication<-as.factor(cols$replication)
  cols$group<-as.factor(paste0(cols$baits,cols$conditions)) 
  if(norm_plots){
    source(paste0(location_folders,"/normalization_plots.R"))
    normalization_plots(sample_names, raw_counts, cols, out_dir)
  }
  
  ### lib size normalization
  dds=DESeqDataSetFromMatrix(countData = raw_counts, colData = cols, design = ~ group)
  sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
  dds <- DESeq(dds, quiet = TRUE, betaPrior=FALSE)
  saveRDS(dds, paste0(output_location,"dds.RDS"))
  
  #run enrich and spec
  source(paste0(location_folders,"/enrichment_score_newInt.R"))
  #dds<- readRDS(paste0(output_location,"dds.RDS"))
  enrichment_score_all_rel<- calculate_enrichment_score(dds, sample_names, p_val_thr, output_location)
  source(paste0(location_folders,"/spec_score_newInt.R"))
  spec_score_all_rel<- calc_spec_score(dds, sample_names, p_val_thr_s, fc_thr, output_location)
  
  #combine
  in_frame_results_sum<- readRDS(paste0(output_location,"in_frame_results_sum_rel_np_freq.RDS"))
  enrichment_score_all_rel<- readRDS(paste0(output_location,"enrichment_score_all_rel.RDS"))
  spec_score_all_rel<- readRDS(paste0(output_location,"spec_score_all_rel.RDS"))
  colnameIF<- "total_max_freq_score_np"
  
  total_scores <- merge(enrichment_score_all_rel[, c("gene","bait","total_bait_enrich_score")], 
                        spec_score_all_rel[,c("gene","bait","total_spec_score")], by=c("gene", "bait"), all = T)
  
  total_scores <- merge(total_scores, in_frame_results_sum[,c("gene","bait",colnameIF)],
                        by=c("gene", "bait"), all = T)
  
  total_scores$sum_scores <- rowSums(total_scores[,c("total_bait_enrich_score", "total_spec_score",
                                                     colnameIF)], na.rm=TRUE)
  saveRDS(total_scores, paste0(output_location,"total_scores.RDS"))
  
}


