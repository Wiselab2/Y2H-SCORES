# get raw reads from Yachie paper lab 
#bait DB
#prey AD
ccc_1_S<- read.table("BFG_Y2H_Yachie/Yachie_Petsalaki_Data_S3/CCC_1_-His.tsv",
                     sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
ccc_1_NS<- read.table("BFG_Y2H_Yachie/Yachie_Petsalaki_Data_S3/CCC_1_+His.tsv",
                     sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)

save_raw_counts<- function(table_s, table_n, bait, prey_list, num_reps= 3, path, mul_counts=T){
  table_s<- table_s[table_s$DB_SYMB==bait,c("AD_SYMB", "COUNT")]
  table_n<- table_n[table_n$DB_SYMB==bait,c("AD_SYMB", "COUNT")]
  table<- data.frame(prey=prey_list, s_col=0, ns_col=0)
  table[match(table_s$AD_SYMB, table$prey), "s_col"]<- table_s$COUNT
  table[match(table_n$AD_SYMB, table$prey), "ns_col"]<- table_n$COUNT
  
  raw_counts<-matrix(c(rep(table[,"s_col"], num_reps), rep(table[,"ns_col"], num_reps)), nrow = length(prey_list), ncol = 2*num_reps)
  raw_counts<- round(raw_counts)
  for(col in 1:ncol(raw_counts)){
    if(mul_counts){
      raw_counts[,col]<- 1000*raw_counts[,col] + round(sample(1:num_reps, 1))
    }
  }
  rownames(raw_counts)<- table$prey
  write.table(raw_counts, paste0(path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)

}

make_input_arguments_file<- function(draft_table, bait_name, num_reps, method_folder="BFG_Y2H_Yachie"){
  draft_table[2,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "S.fastq"), collapse = ";")
  draft_table[4,"argument_value"]<- paste(paste0("/", bait_name, "R", 1:num_reps, "N.fastq"), collapse = ";")
  draft_table[7,"argument_value"]<- paste0(method_folder, "/formated_files_for_scores/", bait_name)
  write.csv(draft_table, paste0("", method_folder, 
                                "/formated_files_for_scores/", bait_name, "_input_arguments.csv"), row.names = F, quote = F)
}

draft_table<- read.csv("draft_input_arguments.csv")
dir_path<- "BFG_Y2H_Yachie/formated_files_for_scores/"
prey_list<- unique(c(ccc_1_NS$AD_SYMB, ccc_1_S$AD_SYMB))
validated_interactions<- read.csv("BFG_Y2H_Yachie/validated_interactions_BFG_Y2H_CCC_extra.csv")
bait_list<- unique(validated_interactions$bait)
fofn<- list()
for(bait in bait_list){
  bait_folder<- paste0("BFG_Y2H_Yachie/formated_files_for_scores/", bait, "/")
  dir.create(bait_folder)
  save_raw_counts(table_s=ccc_1_S, table_n=ccc_1_NS, bait=bait, prey_list, num_reps= 3, path= bait_folder)
  make_input_arguments_file(draft_table, bait_name=bait, num_reps=3, method_folder="BFG_Y2H_Yachie")
  fofn[[bait]]<- paste0("BFG_Y2H_Yachie/formated_files_for_scores/", bait, "_input_arguments.csv")
}
fofn_table<- t(bind_rows(fofn)) 
write.table(fofn_table, "BFG_Y2H_Yachie/formated_files_for_scores/fofn_for_compute_scores_BFG_Y2H.txt", 
            col.names = F, row.names = F, quote = F)

#biogrid
biogrid <- read.table(file="BIOGRID-ALL-4.1.190.tab2.txt", 
                      sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(biogrid)
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",]
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
net<- biogrid[(biogrid$Official.Symbol.Interactor.A %in% bait_list & biogrid$Official.Symbol.Interactor.B %in% prey_list)|
                (biogrid$Official.Symbol.Interactor.A %in% prey_list & biogrid$Official.Symbol.Interactor.B %in% bait_list),]

net<- net[!duplicated(t(apply(net, 1, sort))), ]
net$bait<- ifelse(net$Official.Symbol.Interactor.A %in% bait_list, net$Official.Symbol.Interactor.A, net$Official.Symbol.Interactor.B)
net$prey<- ifelse(net$Official.Symbol.Interactor.A == net$bait, net$Official.Symbol.Interactor.B, net$Official.Symbol.Interactor.A)
net<- net[,c("bait", "prey")]
net$interactor<-1
sum(net$prey %in%  total_scores$prey)
sum(total_scores$prey %in% net$prey)

validated_bfg<- read.csv("BFG_Y2H_Yachie/validated_interactions_BFG_Y2H_CCC_extra.csv")

validated_bfg$interactor<- ifelse(is.na(validated_bfg$Union), validated_bfg$Pairwise_Y2H, validated_bfg$Union)
validated_bfg$interactor[validated_bfg$interactor=="Autoactivator"]<-0
validated_bfg$interactor<- as.numeric(validated_bfg$interactor)

validated_bgf_net<- rbind(validated_bfg[, c(3,4,12)], net)
validated_bgf_net<- validated_bgf_net[!duplicated(validated_bgf_net),]


#RUN SCORES
#Rscript run_scores.R --fofn "BFG_Y2H_Yachie/formated_files_for_scores/fofn_for_compute_scores_BFG_Y2H.txt" --out_dir "BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 1

total_scores<- readRDS("BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/total_scores.RDS")
validated_bgf_net<- readRDS("BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/validated_bfg_biogrid.RDS")

total_scores<- merge(total_scores, validated_bgf_net, by=c("bait","prey"), all.x = T)
total_scores<- merge(total_scores, validated_bfg[,3:5], by=c("bait","prey"), all.x = T)
total_scores[is.na(total_scores$Interaction_Score), "Interaction_Score"]<-0
total_scores[is.na(total_scores$interactor) & total_scores$Interaction_Score>0,"interactor"]<-0

saveRDS(total_scores, "BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/total_scores_validated.RDS")
saveRDS(validated_bgf_net, "BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/validated_bfg_biogrid.RDS")

#plot performance
pdf("BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/ROC_biogrid.pdf", width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
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

basicplot3<- ggplot(total_scores, aes(m = Interaction_Score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot3 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC BGF-Y2H Score"))

total_scores2<- total_scores[!is.na(total_scores$interactor),]
library(PRROC)
pr_2<-pr.curve(scores.class0 = total_scores2$Borda_scores, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_2)

pr_3<-pr.curve(scores.class0 = total_scores2$Specificity_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_3)

pr_1<-pr.curve(scores.class0 = total_scores2$Enrichment_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_1)

pr_4<-pr.curve(scores.class0 = total_scores2[!is.na(total_scores2$Interaction_Score), "Interaction_Score"], 
               weights.class0 = total_scores2[!is.na(total_scores2$Interaction_Score), "interactor"], curve=T)
plot(pr_4)

dev.off()
#put auc in the scenario table
auc<- data.frame(ROC_AUC_Borda=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), 
                 ROC_AUC_Enrichment=round(calc_auc(basicplot)$AUC, 2), ROC_AUC_BGF_Y2H=round(calc_auc(basicplot3)$AUC, 2),
                 PR_AUC_Borda=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, 
                 PR_AUC_Enrichment=pr_1$auc.integral, PR_AUC_BGF_Y2H=pr_4$auc.integral, stringsAsFactors = F)


auc$sample<- "BGF_Y2H"

write.csv(auc, "BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/auc_ccc_biogrid.csv")
write.csv(total_scores, "BFG_Y2H_Yachie/output_BFG_Y2H_CCC_groups_random/total_scores_validated_biogrid.csv")



