# get raw reads from goossens lab 

save_raw_counts<- function(path, save_path){
  raw_counts_se<-data.frame()
  replicates<- list.dirs(path = path, full.names = TRUE, recursive = F)
  for (replicate in replicates){
    counts_per_sample<-read.table(paste0(replicate,"/quant.genes.sf"), row.names = 1, header = T)
    if(nrow(raw_counts_se)==0)
    {
      raw_counts_se<-data.frame(counts_per_sample$NumReads)

    }else{
      raw_counts_se<-cbind(raw_counts_se, counts_per_sample$NumReads)
    }
  }
  raw_counts_se<- round(raw_counts_se)
  rownames(raw_counts_se)<- rownames(counts_per_sample)
  raw_counts_pe<- raw_counts_se[, seq(2,ncol(raw_counts_se), 2)]
  write.table(raw_counts_se, paste0(save_path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)
  write.table(raw_counts_pe, paste0(save_path,strsplit(path, "/")[[1]][length(strsplit(path, "/")[[1]])],"_salmon_counts_PE.matrix"), row.names = T, col.names = F, sep = "\t", quote = F)
  
}

path<- "goossens_lab/raw salmon counts/NINJA"
save<- "goossens_lab/formated_files_for_scores/NINJA/"
save_raw_counts(path, save)

path<- "goossens_lab/raw salmon counts/EMPTY"
save<- "goossens_lab/formated_files_for_scores/EMPTY/"
save_raw_counts(path, save)

#lets add the gene descriptions to the total scores
ann<- read.table("goossens_lab/mart_export.txt",
                  header=T, fill=T, sep = "\t",quote = "\"",stringsAsFactors = F)
ann$prey<-ann[,1]
ann$Description<- paste0(ann$Gene.name, " ", ann$Interpro.Description)
ann<- ann[,c("prey", "Description")]
ann<- ann[!duplicated(ann$prey),]

#RUN SCORES
#Rscript run_scores.R --fofn "goossens_lab/formated_files_for_scores/fofn_for_compute_scores_SE.txt" --out_dir "goossens_lab/output_goossens_SE2/" --spec_p_val 1 --spec_fold_change 0 --enrich_p_val 0 --normalized F

# get previous interactions
library(tidyverse)
taxon_ids <- c("3702")
#biogrid
biogrid <- read.table(file="BIOGRID-ALL-4.1.190.tab2.txt", 
                      sep = '\t',header = TRUE, quote = "\"", fill=TRUE, stringsAsFactors = F)
colnames(biogrid)
biogrid <- biogrid[biogrid$Experimental.System.Type=="physical",c("Organism.Interactor.A", "Organism.Interactor.B",
                                                                  "Systematic.Name.Interactor.A", "Systematic.Name.Interactor.B")]
biogrid<- biogrid[biogrid$Organism.Interactor.A %in% taxon_ids &
                    biogrid$Organism.Interactor.B %in% taxon_ids,c(3,4)]
ninja_net<- biogrid[biogrid$Systematic.Name.Interactor.A=="AT4G28910" |
                      biogrid$Systematic.Name.Interactor.B =="AT4G28910",]
ninja_net<- ninja_net[!duplicated(t(apply(ninja_net, 1, sort))), ]
ninja_net$prey<- apply(ninja_net, 1, FUN=function(x){ifelse(x[1]=="AT4G28910", x[2], x[1])})
ninja_net$bait<- "NINJA"
ninja_net<- ninja_net[,c("bait", "prey")]
#ninja_net<- rbind(ninja_net, data.frame(bait=c("NINJA","NINJA"), prey=c("AT3G17860","AT4G14713")))
ninja_net<- rbind(ninja_net, data.frame(bait="NINJA", prey="AT3G17860"))
ninja_net$interactor<- 1
ninja_net<- rbind(ninja_net, data.frame(bait="NINJA", prey=c("AT1G34340", "AT3G06850", "AT4G36480", "AT5G47810"), interactor=0))

ninja_net_SE<- ninja_net
ninja_net_SE$bait<- "NINJA_SE"
saveRDS(ninja_net_SE, "goossens_lab/output_goossens_SE2/ninja_net.RDS")
ninja_net_SE<- readRDS("goossens_lab/output_goossens_SE2/ninja_net.RDS")

#For SE
total_scores<- readRDS("goossens_lab/output_goossens_SE2/total_scores.RDS")
#total_scores<- merge(total_scores, ann, by="prey", all.x = T)
total_scores<- merge(total_scores, ninja_net_SE, by=c("bait","prey"), all.x = T)
saveRDS(total_scores, "goossens_lab/output_goossens_SE2/total_scores_validated.RDS")

#plot performance
table3_paper<- read.csv("goossens_lab/goosens_table3_paper.csv")
total_scores_SE<- readRDS("goossens_lab/output_goossens_SE2/total_scores_validated.RDS")
total_scores_SE<- total_scores_SE[total_scores_SE$bait=="NINJA_SE",]
total_scores_SE<- merge(total_scores_SE, table3_paper[,c(2,6)], by.x = "prey", by.y="Gene_ID", all.x=T)
total_scores_SE[is.na(total_scores_SE$SNRFPKM), "SNRFPKM"]<- 0
total_scores_SE[is.na(total_scores_SE$interactor) & total_scores_SE$SNRFPKM>0,"interactor"]<-0
total_scores_SE<- total_scores_SE[(total_scores_SE$interactor==1 & !is.na(total_scores_SE$interactor)) |!is.na(total_scores_SE$SNRFPKM),]
total_scores_SE[is.na(total_scores_SE$SNRFPKM) , "SNRFPKM"]<- 0
total_scores_SE[is.na(total_scores_SE)]<- 0

pdf("goossens_lab/output_goossens_SE/ROC2.pdf", width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
library(plotROC)
basicplot1<- ggplot(total_scores_SE, aes(m = Borda_scores, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot1 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot1)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC Borda Score"))

basicplot2<- ggplot(total_scores_SE, aes(m = Specificity_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot2 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot2)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC Specificity Score"))

basicplot<- ggplot(total_scores_SE, aes(m = In_frame_score, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC In-frame Score"))

basicplot3<- ggplot(total_scores_SE, aes(m = SNRFPKM, d = interactor))+ geom_roc(n.cuts=8,labels=F) 
print(basicplot3 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 18) + theme_grey(base_size = 70)+
        annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 22)  + 
        theme(axis.text.x = element_text(angle = 90)) + ggtitle("ROC SNR FPKM Score"))

total_scores2<- total_scores_SE[!is.na(total_scores_SE$interactor),]

library(PRROC)
pr_2<-pr.curve(scores.class0 = total_scores2$Borda_scores, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_2)

pr_3<-pr.curve(scores.class0 = total_scores2$Specificity_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_3)

pr_1<-pr.curve(scores.class0 = total_scores2$In_frame_score, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_1)

pr_4<-pr.curve(scores.class0 = total_scores2$SNRFPKM, weights.class0 = total_scores2$interactor, curve=T)
plot(pr_4)
dev.off()
#put auc in the scenario table
auc_SE<- data.frame(ROC_AUC_Borda=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), 
                 ROC_AUC_in_frame=round(calc_auc(basicplot)$AUC, 2), ROC_AUC_SNR=round(calc_auc(basicplot3)$AUC, 2),
                 PR_AUC_Borda=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, 
                 PR_AUC_in_frame=pr_1$auc.integral, PR_AUC_SNR=pr_4$auc.integral, stringsAsFactors = F)

write.csv(auc_SE, "goossens_lab/output_goossens_SE/auc_SE2.csv")
write.csv(total_scores_SE, "goossens_lab/output_goossens_SE/total_scores_SE2_goosens.csv")
total_scores_SE<- read.csv("goossens_lab/output_goossens_SE/total_scores_SE2_goosens.csv")


