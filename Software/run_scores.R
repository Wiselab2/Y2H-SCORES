require(psych)
require(reshape2)
require(DESeq2)
require(tidyverse)
#require(miRLAB)
require(optparse)
#require(doParallel)
#parameters to be changed by the user


option_list = list(
  make_option("--fofn", type="character", default=NULL, 
              help="Enter a list of csv file names for each bait to be scored.", metavar="character"),
  make_option("--spec_groups", type="character", default=NULL, 
              help="Enter a list of sample subsets that should be tested together for specificity.", metavar="character"),
  make_option("--spec_p_val", type="double", default=1, 
              help="The maximum cut off for p values for specificity score [default= %default]"),
  make_option("--spec_fold_change", type="double", default=0, 
              help="Minimum fold change for specificity score [default= %default]"),
  make_option("--enrich_p_val", type="double", default=1, 
              help="The maximum cut off for p values for enrichment score [default= %default]"),
  make_option("--enrich_fold_change", type="double", default=0, 
              help="Minimum fold change for enrichment score [default= %default]"),
  make_option("--out_dir", type="character", default=NULL, 
              help="Enter the name of the output directory where all the results will be stored."),
  make_option("--normalized", type="logical", default=F, 
              help="Enter if the count files are already normalized [default= %default]. Default implements lib size normalization")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$fofn)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input fofn).n", call.=FALSE)
}

if (is.null(opt$out_dir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input out_dir).n", call.=FALSE)
}


#args <- commandArgs()
#print(args)
# p_val_thr_s<-as.numeric(args[6])
# fc_thr_s<-as.numeric(args[7])
# p_val_thr<-as.numeric(args[8])
# out_dir<-strsplit(args[9],"::")[[1]]
# spec_groups_path<- strsplit(args[10],"::")[[1]]
# salmon_counts_list<-strsplit(args[11],"::")[[1]]
# final_reports_list<-strsplit(args[12],"::")[[1]]
# sample_names<-strsplit(args[13],"::")[[1]]
# s_num_replicates<-as.numeric(strsplit(args[14],"::")[[1]])
# n_num_replicates<-as.numeric(strsplit(args[15],"::")[[1]])
# names_s_replicates<-strsplit(args[16],"::")[[1]]


p_val_thr_s<-as.numeric(opt$spec_p_val)
fc_thr_s<-as.numeric(opt$spec_fold_change)
p_val_thr<-as.numeric(opt$enrich_p_val)
fc_thr<-as.numeric(opt$enrich_fold_change)
out_dir<-opt$out_dir
spec_groups<- opt$spec_groups
normalized<- opt$normalized

config_files <- scan(opt$fofn, what = character())
#config_files <- scan("~/iowa_state/lab/Y2H publication/second_submission/paper_datasets/DEEPN/formated_files_for_scores/fofn_for_compute_scores.txt", what = character())
salmon<- list()
report<- list()
s_n_rep<- list()
n_n_rep<- list()
s_name_rep<- list()

for(file in config_files){
  bait_config<- read.csv(file)
  bait_name<- str_split(bait_config[7,"argument_value"], "/")[[1]]
  bait_name<- bait_name[[length(bait_name)]]
  
  salmon[[bait_name]]<- paste0(bait_config[7,"argument_value"], "/", bait_name, "_salmon_counts.matrix")
  report[[bait_name]]<- paste0(bait_config[7,"argument_value"], "/", bait_name, "_final_report.csv")
  s_n_rep[[bait_name]]<- length(str_split(bait_config[2,"argument_value"], ";")[[1]])
  if(is.na(bait_config[4,"argument_value"]) | bait_config[4,"argument_value"]==""){
    n_n_rep[[bait_name]]<- 0
  }else{
    n_n_rep[[bait_name]]<- length(str_split(bait_config[4,"argument_value"], ";")[[1]])
  }
  s_name_rep[[bait_name]]<- str_split(bait_config[2,"argument_value"], ";")
}

salmon_counts_list<-unlist(salmon)
final_reports_list<-unlist(report)
sample_names<-names(salmon)
s_num_replicates<-unlist(s_n_rep)
n_num_replicates<-unlist(n_n_rep)
names_s_replicates<-unlist(s_name_rep)
names_s_replicates<- gsub('.fastq.*','', gsub('.*\\/','',names_s_replicates))


location_folders<- getwd()
dir.create(out_dir)
print("Checking parameters for running Y2H-SCORES")
print(paste0("Specificity p-value threshold: ", p_val_thr_s))
print(paste0("Specificity fold change threshold: ", fc_thr_s))
print(paste0("Enrichment p-value threshold: ", p_val_thr))
print(paste0("Enrichment fold change threshold: ", fc_thr))
print(paste0("Output directory: ", out_dir))
print("Salmon counts:")
print(salmon_counts_list)
print("Fusion counts:")
print(final_reports_list)
print("Sample names:")
print(sample_names)
print("Number of replicates for selected samples:")
print(s_num_replicates)
print(str(s_num_replicates))
print("Number of replicates for non-selected samples:")
print(n_num_replicates)
print(str(n_num_replicates))
print("Path for the file location:")
print(location_folders)
print("Name of replicates for selected samples:")
print(names_s_replicates)
print("Sample groups for specificity score:")
print(spec_groups)
print("Are the samples normalized?")
print(normalized)

############################################################################
output_location<- out_dir
# Read in the data from the files and create a data frame for the entire data
raw_counts<-data.frame()
for(sample_num in 1:length(sample_names)){
  if(s_num_replicates[sample_num]>0 & n_num_replicates[sample_num]>0){
    column_names<-c(paste0("R",1:s_num_replicates[sample_num], "S"), paste0("R",1:n_num_replicates[sample_num], "N"))
  }
  if(n_num_replicates[sample_num]==0){
    column_names<-paste0("R",1:s_num_replicates[sample_num], "S")
  }
  if(s_num_replicates[sample_num]==1 & n_num_replicates[sample_num]==0){
    return(print("Not enough replicates to run Y2H-SCORES"))
  }
  raw_counts_per_sample<-read.table(salmon_counts_list[sample_num], row.names = 1)
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
if(normalized){
  normFactor <- rep(1, ncol(raw_counts))
  names(normFactor)<- names(totCounts)
}else{
  normFactor <- totCounts/geometric.mean(totCounts[totCounts>0], na.rm = T)
}
normalized_counts<- raw_counts/rep(normFactor, each = (nrow(raw_counts)))
normalized_counts<-round(normalized_counts)
saveRDS(normalized_counts, paste0(output_location,"normalized_counts_salmon.RDS"))
saveRDS(normFactor, paste0(output_location,"normFactor_salmon.RDS"))
saveRDS(raw_counts, paste0(output_location,"raw_counts_salmon.RDS"))

if(!is.null(spec_groups)){
  spec_groups <- scan(spec_groups, what = list(name=character()))
  spec_groups <- sapply(strsplit(spec_groups$name," "), FUN=function(x)strsplit(x,","))
  groups<- spec_groups
}else{
  ngroups<- length(sample_names)%/%10 + ifelse(length(sample_names)%%10<3 & length(sample_names)%/%10>0, 0, 1)
  sampling_set<- sample(sample_names, length(sample_names), replace = F)
  groups<- sapply(1:ngroups, FUN=function(x) if(x!=ngroups) list(sampling_set[(10*(x-1)+1):(10*x)]) else list(sampling_set[(10*(x-1)+1):length(sampling_set)]))
  saveRDS(groups, paste0(output_location,"random_spec_groups.RDS"))
}

# normFactor<- readRDS(paste0(output_location,"normFactor_salmon.RDS"))

#run in frame
if(file.exists(final_reports_list[1])){
  source(paste0(location_folders,"/in_frame_score.R"))
  in_frame_results_sum<- calc_in_frame(normFactor, sample_names, final_reports_list, output_location, names_s_replicates)
  saveRDS(in_frame_results_sum, paste0(output_location, "in_frame_results_sum.RDS"))
  
}

# if(sum(n_num_replicates)>0){
#   source(paste0(location_folders,"/in_frame_score.R"))
#   in_frame_results_sum<- calc_in_frame(normFactor, sample_names, final_reports_list, output_location, names_s_replicates)
# }
# Run DESeq2 for enriched and specificity score
condition<-c("S","N")
cols<-data.frame(baits=rep(sample_names, (n_num_replicates+s_num_replicates)), 
                 conditions=unlist(list(apply(data.frame(s_num_replicates, n_num_replicates), 1, FUN=function(x){rep(condition,c(x[1], x[2]))}))),
                 replication=unlist(list(apply(data.frame(s_num_replicates, n_num_replicates), 1, FUN=function(x){ifelse(x[2]>0,return(c(1:x[1],1:x[2])),return(1:x[1]))}))))
rownames(cols)<-colnames(raw_counts)
cols$conditions<-as.factor(cols$conditions)
cols$replication<-as.factor(cols$replication)
cols$group<-as.factor(paste0(cols$baits,cols$conditions)) 
saveRDS(cols, paste0(output_location,"cols.RDS"))

### lib size normalization
run_DESeq2 = function(dds){tryCatch(DESeq(dds, quiet = TRUE, betaPrior=FALSE), 
                                    error = function(e){ 
                                      print(e)
                                      dds <- estimateDispersionsGeneEst(dds)
                                      dispersions(dds) <- mcols(dds)$dispGeneEst
                                      dds<- nbinomWaldTest(dds)})}

dds_list<- list()
for(g in 1:length(groups)){
  dds_list[[g]]<- DESeqDataSetFromMatrix(countData = raw_counts[,rownames(cols[cols$baits %in% groups[[g]],])], 
                                         colData = cols[cols$baits %in% groups[[g]], ], design = ~ group)
  sizeFactors(dds_list[[g]]) <- normFactor[names(normFactor) %in% rownames(cols[cols$baits %in% groups[[g]],])]
  dds_list[[g]]<- run_DESeq2(dds_list[[g]])
}

# dds=DESeqDataSetFromMatrix(countData = raw_counts, colData = cols, design = ~ group)
# sizeFactors(dds) <- normFactor
# dds<- run_DESeq2(dds)
# run_DESeq2 = function(dds){tryCatch(DESeq(dds, quiet = TRUE, betaPrior=FALSE), 
#                                     error = function(e){ 
#                                       print(e)
#                                       dds <- estimateDispersionsGeneEst(dds)
#                                       dispersions(dds) <- mcols(dds)$dispGeneEst
#                                       dds<- nbinomWaldTest(dds)})}



#Now the contrasts between genotypes and timepoints https://support.bioconductor.org/p/67600/#67612
saveRDS(dds_list, paste0(output_location,"dds.RDS"))

#run enrich and spec
if(sum(n_num_replicates)>0){
  source(paste0(location_folders,"/enrichment_score_groups.R"))
  #dds<- readRDS(paste0(output_location,"dds.RDS"))
  enrichment_score_all_rel<- calculate_enrichment_score(dds_list, sample_names, p_val_thr, fc_thr, output_location, groups)
  saveRDS(enrichment_score_all_rel, paste0(output_location, "enrichment_score_all_rel.RDS"))
}
if(length(salmon_counts_list)>1){
  #source(paste0(location_folders,"/spec_score.R"))
  #  spec_score_all_rel<- calc_spec_score(dds, sample_names, p_val_thr_s, fc_thr_s, output_location)
  source(paste0(location_folders,"/spec_score_groups.R"))
  spec_score_all_rel<- calc_spec_score(dds_list, sample_names, p_val_thr_s, fc_thr_s, output_location, groups)
  saveRDS(spec_score_all_rel, paste0(output_location, "spec_score_all_rel.RDS"))
}

#combine
# in_frame_results_sum<- readRDS(paste0(output_location,"in_frame_results_sum.RDS"))
# enrichment_score_all_rel<- readRDS(paste0(output_location,"enrichment_score_all_rel.RDS"))
# spec_score_all_rel<- readRDS(paste0(output_location,"spec_score_all_rel.RDS"))

Borda<- function (listCEmatrices) 
{
  noMethods = length(listCEmatrices)
  effectnames = rownames(listCEmatrices[[1]])
  causenames = colnames(listCEmatrices[[1]])
  res = matrix(nrow = length(effectnames), ncol = length(causenames), 
               0)
  colnames(res) = causenames
  rownames(res) = effectnames
  for (i in 1:noMethods) {
    for (j in 1:length(causenames)) {
      cormat = listCEmatrices[[i]][, j]
      cormat = cormat[order(-abs(cormat))]
      for (k in 1:length(cormat)) {
        cormat[k] = k
      }
      rn = names(cormat)
      for (gene in rn) {
        res[gene, j] = res[gene, j] + cormat[gene]
      }
    }
  }
  res = res/noMethods
  res = length(effectnames)/res
}

colnameIF<- "total_max_freq_score"
if(sum(n_num_replicates)>0 & length(salmon_counts_list)>1 & file.exists(final_reports_list[1])){
  total_scores <- merge(enrichment_score_all_rel[, c("gene","bait","total_bait_enrich_score")], 
                        spec_score_all_rel[,c("gene","bait","total_spec_score")], by=c("gene", "bait"), all = T)
  total_scores <- merge(total_scores, in_frame_results_sum[,c("gene","bait",colnameIF, "transcript")],
                        by=c("gene", "bait"), all = T)
  total_scores$sum_scores <- rowSums(total_scores[,c("total_bait_enrich_score", "total_spec_score",
                                                     colnameIF)], na.rm=TRUE)
  #borda counts for the scores
  
  total_scores[is.na(total_scores)]<- 0
  rownames(total_scores)<- paste0(total_scores$gene,"_", total_scores$bait)
  enrich_mat<- as.matrix(x=total_scores$total_bait_enrich_score)
  row.names(enrich_mat)<- row.names(total_scores)
  colnames(enrich_mat)<- "conf"
  spec_mat<- as.matrix(x=total_scores$total_spec_score)
  row.names(spec_mat)<- row.names(total_scores)
  colnames(spec_mat)<- "conf"
  inframe_mat<- as.matrix(x=total_scores$total_max_freq_score)
  row.names(inframe_mat)<- row.names(total_scores)
  colnames(inframe_mat)<- "conf"
  borda_scores<- Borda(list(enrich_mat, spec_mat, inframe_mat))
  total_scores$borda<- borda_scores[,1]
  colnames(total_scores)<- c("prey", "bait", "Enrichment_score", "Specificity_score", "In_frame_score",
                             "In_frame_prey_transcripts", "Sum_scores", "Borda_scores")
}
if(sum(n_num_replicates)>0 & length(salmon_counts_list)==1 & file.exists(final_reports_list[1])){
  total_scores <- merge(enrichment_score_all_rel[, c("gene","bait","total_bait_enrich_score")], 
                        in_frame_results_sum[,c("gene","bait",colnameIF, "transcript")], by=c("gene", "bait"), all = T)
  total_scores$sum_scores <- rowSums(total_scores[,c("total_bait_enrich_score", colnameIF)], na.rm=TRUE)
  #borda counts for the scores
  
  total_scores[is.na(total_scores)]<- 0
  rownames(total_scores)<- paste0(total_scores$gene,"_", total_scores$bait)
  enrich_mat<- as.matrix(x=total_scores$total_bait_enrich_score)
  row.names(enrich_mat)<- row.names(total_scores)
  colnames(enrich_mat)<- "conf"
  inframe_mat<- as.matrix(x=total_scores$total_max_freq_score)
  row.names(inframe_mat)<- row.names(total_scores)
  colnames(inframe_mat)<- "conf"
  borda_scores<- Borda(list(enrich_mat, inframe_mat))
  total_scores$borda<- borda_scores[,1]
  
  colnames(total_scores)<- c("prey", "bait", "Enrichment_score", "In_frame_score",
                             "In_frame_prey_transcripts", "Sum_scores", "Borda_scores")
}
if(sum(n_num_replicates)==0 & length(salmon_counts_list)>1 & file.exists(final_reports_list[1])){
  total_scores <- merge(spec_score_all_rel[,c("gene","bait","total_spec_score")], 
                        in_frame_results_sum[,c("gene","bait",colnameIF, "transcript")], by=c("gene", "bait"), all = T)
  total_scores$sum_scores <- rowSums(total_scores[,c("total_spec_score", colnameIF)], na.rm=TRUE)
  #borda counts for the scores
  
  total_scores[is.na(total_scores)]<- 0
  rownames(total_scores)<- paste0(total_scores$gene,"_", total_scores$bait)
  spec_mat<- as.matrix(x=total_scores$total_spec_score)
  row.names(spec_mat)<- row.names(total_scores)
  colnames(spec_mat)<- "conf"
  inframe_mat<- as.matrix(x=total_scores$total_max_freq_score)
  row.names(inframe_mat)<- row.names(total_scores)
  colnames(inframe_mat)<- "conf"
  borda_scores<- Borda(list(spec_mat, inframe_mat))
  total_scores$borda<- borda_scores[,1]
  
  colnames(total_scores)<- c("prey", "bait", "Specificity_score", "In_frame_score",
                             "In_frame_prey_transcripts", "Sum_scores", "Borda_scores")
}
if(sum(n_num_replicates)>0 & length(salmon_counts_list)>1 & !file.exists(final_reports_list[1])){
  total_scores <- merge(enrichment_score_all_rel[, c("gene","bait","total_bait_enrich_score")], 
                        spec_score_all_rel[,c("gene","bait","total_spec_score")], by=c("gene", "bait"), all = T)
  
  total_scores$sum_scores <- rowSums(total_scores[,c("total_bait_enrich_score", "total_spec_score")], na.rm=TRUE)
  #borda counts for the scores
  
  total_scores[is.na(total_scores)]<- 0
  rownames(total_scores)<- paste0(total_scores$gene,"_", total_scores$bait)
  enrich_mat<- as.matrix(x=total_scores$total_bait_enrich_score)
  row.names(enrich_mat)<- row.names(total_scores)
  colnames(enrich_mat)<- "conf"
  spec_mat<- as.matrix(x=total_scores$total_spec_score)
  row.names(spec_mat)<- row.names(total_scores)
  colnames(spec_mat)<- "conf"
  borda_scores<- Borda(list(enrich_mat, spec_mat))
  total_scores$borda<- borda_scores[,1]
  colnames(total_scores)<- c("prey", "bait", "Enrichment_score", "Specificity_score", "Sum_scores", "Borda_scores")
}
write.csv(total_scores, paste0(output_location,"Total_scores.csv"), row.names = F)
saveRDS(total_scores, paste0(output_location,"total_scores.RDS"))

