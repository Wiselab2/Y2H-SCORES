require(psych)
require(reshape2)
require(DESeq2)
#require(doParallel)
#parameters to be changed by the user
args <- commandArgs()
#print(args)
p_val_thr_s<-as.numeric(args[6])
fc_thr<-as.numeric(args[7])
p_val_thr<-as.numeric(args[8])
out_dir<-strsplit(args[9],"::")[[1]]
salmon_counts_list<-strsplit(args[10],"::")[[1]]
final_reports_list<-strsplit(args[11],"::")[[1]]
sample_names<-strsplit(args[12],"::")[[1]]
num_replicates<-as.numeric(strsplit(args[13],"::")[[1]])
location_folders<- getwd()
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
#saveRDS(normalized_counts, paste0(output_location,"normalized_counts_salmon.RDS"))
#saveRDS(normFactor, paste0(output_location,"normFactor_salmon.RDS"))
#saveRDS(raw_counts, paste0(output_location,"raw_counts_salmon.RDS"))
#run in frame
source(paste0(location_folders,"/in_frame_score.R"))
#normFactor<- readRDS(paste0(output_location,"normFactor_salmon.RDS"))
in_frame_results_sum<- calc_in_frame(normFactor, sample_names, final_reports_list, output_location)
# Run DESeq2 for enriched and specificity score
condition<-c("S","N")
cols<-data.frame(baits=rep(sample_names, 2*num_replicates), 
                 conditions=unlist(list(apply(data.frame(num_replicates), 1, FUN=function(x){rep(condition,each=x[1])}))), 
                 replication=unlist(list(apply(data.frame(num_replicates), 1, FUN=function(x){rep(1:x[1],2)}))))
rownames(cols)<-colnames(raw_counts)
cols$conditions<-as.factor(cols$conditions)
cols$replication<-as.factor(cols$replication)
cols$group<-as.factor(paste0(cols$baits,cols$conditions)) 
#saveRDS(cols, paste0(output_location,"cols.RDS"))

### lib size normalization
dds=DESeqDataSetFromMatrix(countData = raw_counts, colData = cols, design = ~ group)
sizeFactors(dds) <- totCounts/geometric.mean(totCounts)
dds <- DESeq(dds, quiet = TRUE, betaPrior=FALSE)
#Now the contrasts between genotypes and timepoints https://support.bioconductor.org/p/67600/#67612
#saveRDS(dds, paste0(output_location,"dds.RDS"))

#run enrich and spec
source(paste0(location_folders,"/enrichment_score.R"))
#dds<- readRDS(paste0(output_location,"dds.RDS"))
enrichment_score_all_rel<- calculate_enrichment_score(dds, sample_names, p_val_thr, output_location)
source(paste0(location_folders,"/spec_score.R"))
spec_score_all_rel<- calc_spec_score(dds, sample_names, p_val_thr_s, fc_thr, output_location)

#combine
# in_frame_results_sum<- readRDS(paste0(output_location,"in_frame_results_sum.RDS"))
# enrichment_score_all_rel<- readRDS(paste0(output_location,"enrichment_score_all_rel.RDS"))
# spec_score_all_rel<- readRDS(paste0(output_location,"spec_score_all_rel.RDS"))

colnameIF<- "total_max_freq_score"

total_scores <- merge(enrichment_score_all_rel[, c("gene","bait","total_bait_enrich_score")], 
                      spec_score_all_rel[,c("gene","bait","total_spec_score")], by=c("gene", "bait"), all = T)

total_scores <- merge(total_scores, in_frame_results_sum[,c("gene","bait",colnameIF, "transcript")],
                      by=c("gene", "bait"), all = T)

total_scores$sum_scores <- rowSums(total_scores[,c("total_bait_enrich_score", "total_spec_score",
                                                   colnameIF)], na.rm=TRUE)
#borda counts for the scores
library(miRLAB)
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
write.csv(total_scores, paste0(output_location,"total_scores.csv"), row.names = F)
#saveRDS(total_scores, paste0(output_location,"total_scores.RDS"))

