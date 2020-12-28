# get raw reads from deepn lab 
#description of the S3 table from the deepn paper
# 'Base' refers to average ppm value in NON-selected conditions.  This is the average of the non-selected diploid libraries for 2 independent Vector alone populations, the non-selected Rab5-GTP bait population and the Rab5-GDP non-selected population.  
# Value is the number of mapped sequences that were found assigned to each gene per million mapped reads determined (PPM). 
# 'Vec' is the PPM value for the Vector-alone population under selection in SD-His media.  
# Value represents average of 2 independent Vector-alone populations. 
# 'Rab5GDP' is the PPM value for the Rab5 (SN)-bait plasmid population under selection in SD-HIS media.  
# 'Rab5GTP' is the PPM value for the Rab5 (QL)-bait plasmid population under selection in SD-HIS media.  
# 
# Enr(GDP) shows log('Base' PPM/'Rab5GDP' PPM)
# Enr(GTP) shows log('Base' PPM/'Rab5GDP' PPM)
# 
# p(GDP > Vec) is the ranking for whether the given gene is significantly enriched under selection on Rab5GDP bait vs Vector-alone bait
# p(GTP > Vec) is the ranking for whether the given gene is significantly enriched under selection on Rab5GTP bait vs Vector-alone bait
# p(GDP > GTP) is the ranking for whether the given gene is significantly enriched under selection on Rab5GDP bait vs  Rab5GTP bait
# p(GTP > Vec) is the ranking for whether the given gene is significantly enriched under selection on Rab5GTP bait vs Rab5GDP bait
# pGDP is the ranking for how well a given gene enriches on Rab5GDP vs both the Vector alone and the Rab5GTP baits.
# pGTP is the ranking for how well a given gene enriches on Rab5GDP vs both the Vector alone and the Rab5GDP baits.
# Rab5-QL Rab5GTP
# Rab5-SN Rab5GDP

S3_table_deepn_rab<- read.csv("DEEPN/S3_table_deepn_rab.csv")
S3_table_deepn_rab$prey[duplicated(S3_table_deepn_rab$prey)]
S3_table_deepn_rab$prey[duplicated(S3_table_deepn_rab$prey)]<- paste0(S3_table_deepn_rab$prey[duplicated(S3_table_deepn_rab$prey)], "_2")

S3_table_deepn_Ub<- read.csv("DEEPN/S3_table_deepn_Ub.csv")
sum(S3_table_deepn_Ub$prey %in% S3_table_deepn_rab$prey)
S3_table_deepn_Ub$prey[duplicated(S3_table_deepn_Ub$prey)]
S3_table_deepn_Ub$prey[duplicated(S3_table_deepn_Ub$prey)]<- paste0(S3_table_deepn_Ub$prey[duplicated(S3_table_deepn_Ub$prey)], "_2")

S3_table_deepn_Ub4<- read.csv("DEEPN/S3_table_deepn_Ub4.csv")
sum(S3_table_deepn_Ub4$prey %in% S3_table_deepn_rab$prey)
sum(S3_table_deepn_Ub4$prey %in% S3_table_deepn_Ub$prey)

S3_table_deepn_Ub4$prey[duplicated(S3_table_deepn_Ub4$prey)]
S3_table_deepn_Ub4$prey[duplicated(S3_table_deepn_Ub4$prey)]<- paste0(S3_table_deepn_Ub4$prey[duplicated(S3_table_deepn_Ub4$prey)], "_2")

S3_table_deepn<- merge(S3_table_deepn_rab, S3_table_deepn_Ub, by="prey", all=T)
S3_table_deepn<- merge(S3_table_deepn, S3_table_deepn_Ub4, by="prey", all=T)
S3_table_deepn[is.na(S3_table_deepn)]<- 0
write.csv(S3_table_deepn, "DEEPN/S3_table_deepn.csv", row.names = F)

save_raw_counts<- function(table, ns_col, s_col, num_reps= 3, save_path){
  raw_counts<-matrix(c(rep(table[,s_col], num_reps), rep(table[,ns_col], num_reps)), nrow = nrow(table), ncol = 2*num_reps)
  raw_counts<- round(raw_counts)
  for(col in 2:ncol(raw_counts)){
    raw_counts[,col]<- raw_counts[,col] + round(sample(1:ncol(raw_counts), 1))
  }
  rownames(raw_counts)<- table$prey
  write.table(raw_counts, paste0(save_path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)

}
path<- "DEEPN/formated_files_for_scores/EMPTY/"
save_raw_counts(table=S3_table_deepn, ns_col="Base", s_col= "Vec", num_reps= 3, save_path=path)

path<- "DEEPN/formated_files_for_scores/Rab5GDP/"
save_raw_counts(table=S3_table_deepn, ns_col="Base", s_col= "Rab5GDP", num_reps= 3, save_path=path)

path<- "DEEPN/formated_files_for_scores/Rab5GTP/"
save_raw_counts(table=S3_table_deepn, ns_col="Base", s_col= "Rab5GTP", num_reps= 3, save_path=path)

path<- "DEEPN/formated_files_for_scores/diUb1/"
save_raw_counts(table=S3_table_deepn, ns_col="Base_diUb", s_col= "diUb1", num_reps= 3, save_path=path)

path<- "DEEPN/formated_files_for_scores/di4spUb/"
save_raw_counts(table=S3_table_deepn, ns_col="Base_diUb4", s_col= "di4spUb", num_reps= 3, save_path=path)

#lets add the gene descriptions to the total scores
# get previous interactions
library(tidyverse)
#taxon_ids <- c("10090")
#biogrid
biogrid <- read.table(file="BIOGRID-ALL-4.1.190.tab2.txt", 
                      sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(biogrid)
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",]
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")]
#baits<- c("RAB5", "UBI")
rab5_net<- biogrid[unique(c(grep("RAB5", biogrid$Official.Symbol.Interactor.A),
                      grep("RAB5", biogrid$Official.Symbol.Interactor.B))),]
rab5_net<- rab5_net[!duplicated(t(apply(rab5_net, 1, sort))), ]
rab5_net$prey<- apply(rab5_net, 1, FUN=function(x){ifelse(grepl("RAB5", x[1]), paste(toupper(substr(x[2], 1, 1)), tolower(substr(x[2], 2, nchar(x[2]))), sep=""),
                                                          paste(toupper(substr(x[1], 1, 1)), tolower(substr(x[1], 2, nchar(x[1]))), sep=""))})
rab5_net$bait<- "Rab5GTP"
rab5_net<- rab5_net[,c("bait", "prey")]
rab5_net$interactor<- 1

ub_net<- biogrid[unique(c(grep("UB", biogrid$Official.Symbol.Interactor.A),
                            grep("UB", biogrid$Official.Symbol.Interactor.B))),]
ub_net<- ub_net[!duplicated(t(apply(ub_net, 1, sort))), ]
ub_net$prey<- apply(ub_net, 1, FUN=function(x){ifelse(grepl("UB", x[1]), paste(toupper(substr(x[2], 1, 1)), tolower(substr(x[2], 2, nchar(x[2]))), sep=""),
                                                          paste(toupper(substr(x[1], 1, 1)), tolower(substr(x[1], 2, nchar(x[1]))), sep=""))})
ub_net$bait<- "di4spUb"
ub_net<- ub_net[,c("bait", "prey")]
ub_net$interactor<- 1
write.csv(ub_net, "DEEPN/output_DEEPN_groups_norm/ub_net.csv")
write.csv(rab5_net, "DEEPN/output_DEEPN_groups_norm/rab5_net.csv")
rab5_net<- read.csv("DEEPN/output_DEEPN_groups_norm/rab5_net.csv", row.names = 1)
ub_net<- read.csv("DEEPN/output_DEEPN_groups_norm/ub_net.csv", row.names = 1)

sum(unique(rab5_net$prey) %in% S3_table_deepn$prey)
sum(unique(ub_net$prey) %in% S3_table_deepn$prey)

#RUN SCORES
#Rscript run_scores.R --fofn "DEEPN/formated_files_for_scores/fofn_for_compute_scores.txt" --out_dir "DEEPN/output_DEEPN_groups_norm_all/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 1 --normalized T --spec_groups "DEEPN/formated_files_for_scores/spec_groups_deepn.txt"

total_scores<- readRDS("DEEPN/output_DEEPN_groups_norm_all/total_scores.RDS")
total_scores<- total_scores[total_scores$bait %in% c("Rab5GTP", "di4spUb"),]
S3_table_deepn<- read.csv("DEEPN/S3_table_deepn.csv")
scores_deepn<- S3_table_deepn[,c(1,14,15,21,27)]

colnames(scores_deepn)<- c("prey", "Rab5GDP", "Rab5GTP", "diUb1", "di4spUb")
library(reshape2)
scores_deepn<- melt(scores_deepn, variable.name = "bait", value.name = "DEEPN_score")
total_scores<- merge(total_scores, scores_deepn, by=c("bait","prey"), all.x = T)

validated_deepn<- read.csv("DEEPN/validated_interactions_deepn.csv")
validated_deepn_biogrid<- rbind(validated_deepn, rab5_net, ub_net)
library(tidyverse)
validated_deepn_biogrid<- validated_deepn_biogrid %>% group_by(bait,prey) %>% summarise(interactor=min(interactor))
validated_deepn_biogrid<- validated_deepn_biogrid[!duplicated(validated_deepn_biogrid),]
total_scores<- merge(total_scores, validated_deepn_biogrid, by=c("bait","prey"), all.x = T)
total_scores[is.na(total_scores$interactor) & total_scores$DEEPN_score>0,"interactor"]<-0

saveRDS(total_scores, "DEEPN/output_DEEPN_groups_norm_all/total_scores_validated.RDS")
saveRDS(validated_deepn_biogrid, "DEEPN/output_DEEPN_groups_norm_all/validated_deepn_biogrid2.RDS")


total_scores_all<- readRDS("DEEPN/output_DEEPN_groups_norm_all/total_scores_validated.RDS")
rab5<- unique(rab5_net$prey)[unique(rab5_net$prey) %in% S3_table_deepn[S3_table_deepn$Base>1 , "prey"]]
ub<- unique(ub_net$prey)[unique(ub_net$prey) %in% S3_table_deepn[S3_table_deepn$Base_diUb4>1, "prey"]]

rab5_scores<- unique(total_scores_all[total_scores_all$bait %in% c("Rab5GTP") & total_scores_all$prey %in% rab5, "prey"])
ub_scores<- unique(total_scores_all[total_scores_all$bait %in% c("di4spUb") & total_scores_all$prey %in% ub, "prey"])


pdf("DEEPN/output_DEEPN_groups_norm_all/ROC.pdf", width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
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

basicplot3<- ggplot(total_scores, aes(m = DEEPN_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot3 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC DEEPN Score"))

total_scores2<- total_scores[!is.na(total_scores$interactor),]
library(PRROC)
pr_2<-pr.curve(scores.class0 = total_scores2$Borda_scores, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_2)

pr_3<-pr.curve(scores.class0 = total_scores2$Specificity_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_3)

pr_1<-pr.curve(scores.class0 = total_scores2$Enrichment_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_1)

pr_4<-pr.curve(scores.class0 = total_scores2[!is.na(total_scores2$DEEPN_score), "DEEPN_score"], 
               weights.class0 = total_scores2[!is.na(total_scores2$DEEPN_score), "interactor"], curve=T)
plot(pr_4)

dev.off()
#put auc in the scenario table
auc<- data.frame(ROC_AUC_Borda=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), 
                 ROC_AUC_Enrichment=round(calc_auc(basicplot)$AUC, 2), ROC_AUC_DEEPN=round(calc_auc(basicplot3)$AUC, 2),
                 PR_AUC_Borda=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, 
                 PR_AUC_Enrichment=pr_1$auc.integral, PR_AUC_DEEPN=pr_4$auc.integral, stringsAsFactors = F)


auc$sample<- "DEEPN"

write.csv(auc, "DEEPN/output_DEEPN_groups_norm_all/auc.csv")
write.csv(total_scores, "DEEPN/output_DEEPN_groups_norm_all/total_scores_DEEPN.csv")


