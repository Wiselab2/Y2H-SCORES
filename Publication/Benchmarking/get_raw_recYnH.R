# get raw reads from rec-YnH paper lab 
#NS -Trp RS
#S many reportes RIS
#use corrected counts
#run scores with lib size norm
#DB bait, AD prey
#EXP are replicates
#validation set in S3 and S4
#interaction scores inavIS
#use 25k files


save_raw_counts<- function(s_reps, n_reps, bait, prey_list, path, file_path){
  s_rep1<- read.table(paste0(file_path, s_reps[1]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  s_rep2<- read.table(paste0(file_path, s_reps[2]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  s_rep3<- read.table(paste0(file_path, s_reps[3]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  
  n_rep1<- read.table(paste0(file_path, n_reps[1]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  n_rep2<- read.table(paste0(file_path, n_reps[2]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  n_rep3<- read.table(paste0(file_path, n_reps[3]), sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)[,bait]
  
  raw_counts<-matrix(c(s_rep1,s_rep2, s_rep3, n_rep1, n_rep2, n_rep3), nrow = length(n_rep1), ncol = 6)
  #raw_counts<- round(raw_counts)
  for(col in 1:ncol(raw_counts)){
    raw_counts[,col]<- 1e5*raw_counts[,col] + round(sample(1:3, 1))
  }
  raw_counts<- round(raw_counts)
  rownames(raw_counts)<- prey_list
  
  write.table(raw_counts, paste0(path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)
}

make_input_arguments_file<- function(draft_table, bait_name, num_reps, method_folder="rec-YnH"){
  draft_table[2,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "S.fastq"), collapse = ";")
  draft_table[4,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "N.fastq"), collapse = ";")
  draft_table[7,"argument_value"]<- paste0(method_folder, "/formated_files_for_scores/", bait_name)
  write.csv(draft_table, paste0(method_folder, "/formated_files_for_scores/", bait_name, "_input_arguments.csv"), row.names = F, quote = F)
}

draft_table<- read.csv("draft_input_arguments.csv")
dir_path<- "rec-YnH/formated_files_for_scores/"

bait_list<- unique(colnames(rec_y2h_table)[-1])
prey_list<- rec_y2h_table$DB.Read.1....AD.Read.2.
s_reps<- c("EXP1_Q.25k", "EXP2_Q.25k", "EXP3_Q.25k")
n_reps<- c("EXP1_W.25k", "EXP2_W.25k", "EXP3_W.25k")
fofn<- list()
for(bait in bait_list){
  bait_folder<- paste0(dir_path, bait, "/")
  dir.create(bait_folder)
  save_raw_counts(s_reps=s_reps, n_reps=n_reps, bait, prey_list=prey_list, path= bait_folder, file_path="rec-YnH/recYnH-master/example/output/")
  make_input_arguments_file(draft_table, bait_name=bait, num_reps=3, method_folder="rec-YnH")
  fofn[[bait]]<- paste0("rec-YnH/formated_files_for_scores/", bait, "_input_arguments.csv")
}
fofn_table<- data.frame(unlist(fofn))
write.table(fofn_table, "rec-YnH/formated_files_for_scores/fofn_for_compute_scores_recYnH.txt", 
            col.names = F, row.names = F, quote = F)


#lets add the gene descriptions to the total scores
#biogrid
biogrid <- read.table(file="BIOGRID-ALL-4.1.190.tab2.txt", sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(biogrid)
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",]
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
#baits<- c("RAB5", "UBI")
net<- biogrid[unique(c(grep(paste(bait_list, collapse="|"), biogrid$Official.Symbol.Interactor.A),
                            grep(paste(bait_list, collapse="|"), biogrid$Official.Symbol.Interactor.B))),]
net<- net[net$Official.Symbol.Interactor.A %in% bait_list & net$Official.Symbol.Interactor.B %in% bait_list, ]
net<- net[!duplicated(t(apply(net, 1, sort))), ]
colnames(net)<- c("bait", "prey")
net$interactor<-1
sum(net$prey %in%  total_scores$prey)
sum(total_scores$prey %in% net$prey)
net$interactor<- 1

#total_scores[is.na(total_scores$interactor),"interactor"]<-0
validated_interactions<- read.csv("rec-YnH/validated_recYnH_tests.csv")[,1:3]
validated_interactions2<- read.csv("rec-YnH/validated_recYnH_databases.csv")
validated_interactions<- rbind(validated_interactions, validated_interactions2)
validated_interactions<- validated_interactions[!duplicated(validated_interactions),]
validated_interactions<- rbind(validated_interactions, net)
validated_interactions<- data.frame(bait= c(validated_interactions$bait, validated_interactions$prey),
                                    prey=c(validated_interactions$prey, validated_interactions$bait),
                                    interactor=c(validated_interactions$interactor, validated_interactions$interactor))

library(tidyverse)

validated_interactions<- validated_interactions %>% group_by(bait, prey) %>% summarise(interactor=max(interactor))
validated_interactions<- read.csv("rec-YnH/output_recYnH_groups_random2/validated_interactions.csv", row.names = 1)

#RUN SCORES
#Rscript run_scores_args.R --fofn "rec-YnH/formated_files_for_scores/fofn_for_compute_scores_recYnH.txt" --out_dir "rec-YnH/output_recYnH_groups_random2/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 1 --enrich_fold_change 0 --normalized F

library(reshape2)
total_scores<- readRDS("rec-YnH/output_recYnH_groups_random/total_scores.RDS")

total_scores<- merge(total_scores, validated_interactions, by=c("bait","prey"), all.x = T)

interaction_score<- read.table("rec-YnH/recYnH-master/example/output/EXP.25k.avgIS",
                           sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)

colnames(interaction_score)[1]<- "prey"
interaction_score<- melt(interaction_score, id.vars = "prey", measure.vars = 2:ncol(interaction_score), variable.name = "bait", value.name = "Interaction_score")

total_scores<- merge(total_scores, interaction_score, by=c("bait","prey"), all.x = T)
total_scores[is.na(total_scores$Interaction_score), "Interaction_score"]<- 0
total_scores[is.na(total_scores$interactor) & total_scores$Interaction_score>0,"interactor"]<-0

#plot performance
pdf("rec-YnH/output_recYnH_groups_random/ROC.pdf", width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
library(plotROC)
basicplot1<- ggplot(total_scores, aes(m = Borda_scores, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot1 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot1)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC Borda Score"))

basicplot2<- ggplot(total_scores, aes(m = Specificity_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot2 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot2)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC Specificity Score"))

basicplot<- ggplot(total_scores, aes(m = Enrichment_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC Enrichment Score"))

basicplot3<- ggplot(total_scores, aes(m = Interaction_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot3 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC rec-YnH Interaction Score"))

total_scores2<- total_scores[!is.na(total_scores$interactor),]
#total_scores2<- total_scores
#total_scores2[is.na(total_scores2$interactor),"interactor"]<-0
library(PRROC)
pr_2<-pr.curve(scores.class0 = total_scores2$Borda_scores, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_2)

pr_3<-pr.curve(scores.class0 = total_scores2$Specificity_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_3)

pr_1<-pr.curve(scores.class0 = total_scores2$Enrichment_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_1)

pr_4<-pr.curve(scores.class0 = total_scores2[!is.na(total_scores2$Interaction_score), "Interaction_score"], 
               weights.class0 = total_scores2[!is.na(total_scores2$Interaction_score), "interactor"], curve=T)
plot(pr_4)

dev.off()
#put auc in the scenario table
auc<- data.frame(ROC_AUC_Borda=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), 
                 ROC_AUC_Enrichment=round(calc_auc(basicplot)$AUC, 2), ROC_AUC_rec_Y2H=round(calc_auc(basicplot3)$AUC, 2),
                 PR_AUC_Borda=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, 
                 PR_AUC_Enrichment=pr_1$auc.integral, PR_AUC_rec_Y2H=pr_4$auc.integral, stringsAsFactors = F)


auc$sample<- "rec-YnH"

write.csv(auc, "rec-YnH/output_recYnH_groups_random/auc.csv")
write.csv(total_scores, "rec-YnH/output_recYnH_groups_random/total_scores_recYnH.csv")
saveRDS(total_scores, "rec-YnH/output_recYnH_groups_random/total_scores_validated.RDS")
write.csv(validated_interactions, "rec-YnH/output_recYnH_groups_random/validated_interactions.csv")

