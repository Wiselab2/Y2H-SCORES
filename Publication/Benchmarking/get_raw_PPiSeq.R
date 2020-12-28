# get raw reads from PPiSeq paper 
#bait will be each treatment
#100 bait_prey will be the preys in tables
#25 reps x3 means 75 reps 
#NS MTX(-) 12th gen or MTX(+), same results as 0 gen
#S MTX(+) at 12th generation
#use corrected counts
#run scores with lib size norm

ppi_table<- read.table("PPiSeq/PPiSeq.tsv",
                     sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(ppi_table)
ppi_table$bait_prey<- paste0(ppi_table$bait, "_", ppi_table$prey)

save_raw_counts<- function(ppi_table, bait, path, time=F){
  if(time){
    raw_counts<- data.frame(rep(ppi_table[ppi_table$condition==bait &ppi_table$methotrexate_selection=="yes","bait_prey"], 2),
                            c(ppi_table[ppi_table$condition==bait &ppi_table$methotrexate_selection=="yes","reads_12_generations_corrected"],
                              ppi_table[ppi_table$condition==bait &ppi_table$methotrexate_selection=="yes","reads_0_generations_corrected"]))

  }else{
    raw_counts<- ppi_table[ppi_table$condition==bait,c("bait_prey", "methotrexate_selection", "reads_12_generations_corrected")]
    
  }
  raw_counts$rep<- NA
  
  for(x in unique(raw_counts$bait_prey)){
    raw_counts[raw_counts$bait_prey==x, "rep"]<- 1:150
  }
  raw_counts<- acast(raw_counts, bait_prey ~ rep, value.var = "reads_12_generations_corrected")
  raw_counts<- round(raw_counts)
  
  write.table(raw_counts, paste0(path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)

}

make_input_arguments_file<- function(draft_table, bait_name, num_reps, method_folder="PPiSeq"){
  draft_table[2,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "S.fastq"), collapse = ";")
  draft_table[4,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "N.fastq"), collapse = ";")
  draft_table[7,"argument_value"]<- paste0(method_folder, "/formated_files_for_scores/", bait_name)
  write.csv(draft_table, paste0("", method_folder, 
                                "/formated_files_for_scores/", bait_name, "_input_arguments.csv"), row.names = F, quote = F)
}

draft_table<- read.csv("draft_input_arguments.csv")
dir_path<- "PPiSeq/formated_files_for_scores/"

bait_list<- unique(ppi_table$condition)
fofn<- list()
for(bait in bait_list){
  bait_folder<- paste0("PPiSeq/formated_files_for_scores/", bait, "/")
  dir.create(bait_folder)
  save_raw_counts(ppi_table, bait, path= bait_folder)
  make_input_arguments_file(draft_table, bait_name=bait, num_reps=75, method_folder="PPiSeq")
  fofn[[bait]]<- paste0("PPiSeq/formated_files_for_scores/", bait, "_input_arguments.csv")
}
fofn_table<- data.frame(unlist(fofn))
write.table(fofn_table, "PPiSeq/formated_files_for_scores/fofn_for_compute_scores_PPiSeq.txt", 
            col.names = F, row.names = F, quote = F)
#biogrid
prots<- c("IMD3", "HOM3", "SHR3", "PRS3", "DBP2", "TPO1", "DST1", "FTR1", "FMP45", "HXT1", "PDR5", "RPB9", "SNQ2") 
biogrid <- read.table(file="BIOGRID-ALL-4.1.190.tab2.txt", 
                      sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(biogrid)
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",]
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
net<- biogrid[unique(c(grep(paste(prots, collapse="|"), biogrid$Official.Symbol.Interactor.A),
                            grep(paste(prots, collapse="|"), biogrid$Official.Symbol.Interactor.B))),]
net<- net[net$Official.Symbol.Interactor.A %in% prots & net$Official.Symbol.Interactor.B %in% prots, ]

net<- net[!duplicated(t(apply(net, 1, sort))), ]
net$prey<- apply(net, 1, FUN=function(x){ifelse(x[1] %in% c("IMD3", "HOM3", "SHR3", "PRS3", "DBP2", "TPO1", "DST1", "FTR1", "FMP45"), 
                                                paste0(x[1], "-F[1,2]_", x[2], "-F[3]"),paste0(x[2], "-F[1,2]_", x[1], "-F[3"))})
sum(net$prey %in%  total_scores$prey)
sum(total_scores$prey %in% net$prey)

net$interactor<- 1
validated_ppi<- read.table("PPiSeq/PPsSeq_fitness.tsv", 
                           sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
validated_ppi$prey<- paste0(validated_ppi$bait, "_", validated_ppi$prey)
validated_ppi[validated_ppi$prey %in% net$prey, "Biogrid"]<- 1
validated_ppi$interactor<- ifelse(rowSums(validated_ppi[,7:11])>0, 1,0)
sum(validated_ppi$prey %in% net$prey)
validated_ppi_net<- rbind(validated_ppi[,c(2,12)], net[,3:4])

#RUN SCORES
#Rscript run_scores.R --fofn "PPiSeq/formated_files_for_scores/fofn_for_compute_scores_PPiSeq.txt" --out_dir "PPiSeq/output_PPiSeq_groups_random/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 1 --enrich_fold_change -50 --normalized F

library(reshape2)
total_scores<- readRDS("PPiSeq/output_PPiSeq_groups_random/total_scores.RDS")
total_scores<- merge(total_scores, validated_ppi_net, by="prey", all.x = T)

fitness_score<- read.table("PPiSeq/PPiSeq_fitness_condition.tsv", 
                           sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
fitness_score$prey<- paste0(fitness_score$bait..F.1.2.., "-F[1,2]_", fitness_score$prey..F.3.., "-F[3]")
fitness_score<- melt(fitness_score[3:8], id.vars = "prey", variable.name = "bait", value.name = "Fitness_Score")
total_scores<- merge(total_scores, fitness_score, by=c("bait","prey"), all.x = T)
total_scores[total_scores$prey=="HOM3-F[1,2]_FPR1-F[3]" & total_scores$bait=="FK506", "interactor"]<-0
total_scores[is.na(total_scores$Fitness_Score), "Fitness_Score"]<- 0

#plot performance
pdf("PPiSeq/output_PPiSeq_groups_random/ROC.pdf", width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
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

basicplot3<- ggplot(total_scores, aes(m = Fitness_Score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot3 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC PPiSEq Fitness Score"))

total_scores2<- total_scores[!is.na(total_scores$interactor),]

library(PRROC)
pr_2<-pr.curve(scores.class0 = total_scores2$Borda_scores, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_2)

pr_3<-pr.curve(scores.class0 = total_scores2$Specificity_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_3)

pr_1<-pr.curve(scores.class0 = total_scores2$Enrichment_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_1)

pr_4<-pr.curve(scores.class0 = total_scores2[!is.na(total_scores2$Fitness_Score), "Fitness_Score"], 
               weights.class0 = total_scores2[!is.na(total_scores2$Fitness_Score), "interactor"], curve=T)
plot(pr_4)

dev.off()
#put auc in the scenario table
auc<- data.frame(ROC_AUC_Borda=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), 
                 ROC_AUC_Enrichment=round(calc_auc(basicplot)$AUC, 2), ROC_AUC_BGF_Y2H=round(calc_auc(basicplot3)$AUC, 2),
                 PR_AUC_Borda=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, 
                 PR_AUC_Enrichment=pr_1$auc.integral, PR_AUC_BGF_Y2H=pr_4$auc.integral, stringsAsFactors = F)


auc$sample<- "PPiSeq"

write.csv(auc, "PPiSeq/output_PPiSeq_groups_random/auc.csv")
write.csv(total_scores, "PPiSeq/output_PPiSeq_groups_random_t0/total_scores_PPiSeq.csv")
saveRDS(total_scores, "PPiSeq/output_PPiSeq_groups_random_t0/total_scores_validated.RDS")


