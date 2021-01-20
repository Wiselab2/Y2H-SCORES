#start simulation

#function to sample the non-interactors
simulate_non_interactors<-function(Q_S_table, simulated_prey_names, num_baits){
  library(data.table)
  observed_preys <- Q_S_table
  n_preys<- length(simulated_prey_names)
  #simulated_prey_names <- prey_names
  #simulated_prey_names <- paste0('prey_', 1:n_preys)
  
  observed_preys <- as.data.table(observed_preys)
  n_sim_baits <- num_baits
  preys <- unique(observed_preys$prey)
  n_real_baits <- length(unique(observed_preys$bait))
  sim_preys <- sample(seq(length(preys)), n_preys, replace = TRUE)
  simmed_preys <- data.table(prey = character(),bait = character(),
                             e_ik = numeric(),phi_ikS = numeric())
  for(i in seq(n_preys)){
    prey_name <- preys[sim_preys[i]]
    prey <- observed_preys[prey == prey_name]
    n_individual_prey_baits <- nrow(prey)
    n <- round((n_sim_baits*n_individual_prey_baits)/n_real_baits)
    index <- replicate(sample(seq(n_individual_prey_baits), n_individual_prey_baits), 
                       n = ceiling(n_sim_baits/n_individual_prey_baits))[seq(n)]
    simmed_preys <- rbindlist(list(simmed_preys,
                                   data.table(prey = simulated_prey_names[i], 
                                              bait = paste0('bait', sample(seq(n_sim_baits))),
                                              e_ik = c(prey$e_ik[index], rep(0, n_sim_baits - n)),
                                              phi_ikS = c(prey$phi_ikS[index], rep(0, n_sim_baits - n)))))
  }
  return(data.frame(simmed_preys))
}

galton_walton_sim_table<- function(QS_table, QN_table, t_vector, L_vector, num_preys, num_rep, M_0=3.84e9, 
                                   ns_gen=4, low_conc=F, true_int, high_phi=F, grow_sat=F, overdispersed=F){
  x<- matrix(nrow=num_preys, ncol=length(L_vector), 
             dimnames = list(unique(QS_table$prey), names(L_vector)))
  z<- matrix(nrow=num_preys, ncol=length(L_vector),
             dimnames = list(unique(QS_table$prey), names(L_vector)))
  count_list<- list()
  for(bait in unique(QS_table$bait)){
    qk<- QN_table[,c("prey","q_k")]
    
    if(low_conc){qk[qk$prey %in% true_int[true_int$bait==bait, "prey"], "q_k"]<- min(qk$q_k)}
    x[, grep(paste0(bait,"R"), colnames(x))]<- t(mapply(function(qk)rbinom(n=length(grep(paste0(bait,"R"), colnames(x))), size=M_0, prob=qk), qk$q_k))
    ns_samples<- paste0(rep(bait, each=num_rep),"R",1:num_rep, "N")
    x[, ns_samples]<- x[, ns_samples]*2^ns_gen
    #z[, ns_samples]<- t(mapply(function(m,phi)rnbinom(num_rep, 1/phi, mu=m), x[,ns_samples]%*%diag(L_vector[ns_samples]),QN_table$phi_kN))
    if(overdispersed){QN_table$phi_kN<- sample(QN_table$phi_kN[QN_table$phi_kN >= quantile(QN_table$phi_kN, 0.9)], length(QN_table$phi_kN), replace = T)}
    z[, ns_samples]<- t(mapply(function(n)ceiling(rnbinom(num_rep, 1/QN_table$phi_kN[n], mu=(x[,ns_samples]%*%diag(L_vector[ns_samples]))[n,])), 1:num_preys))
    
    s_samples<- paste0(rep(bait, each=num_rep),"R",1:num_rep, "S")
    qs<- QS_table[QS_table$bait==bait,]
    if(high_phi){qs[qs$prey %in% true_int[true_int$bait==bait, "prey"], "phi_ikS"]<- sample(qs$phi_ikS[qs$phi_ikS >= quantile(qs$phi_ikS, 0.9)], length(true_int[true_int$bait==bait, "prey"]), replace = T)}
    if(overdispersed){qs$phi_ikS<- sample(qs$phi_ik[qs$phi_ikS >= quantile(qs$phi_ikS, 0.9)], length(qs$phi_ikS), replace = T)}
    
    for(sample in s_samples){
      #bait<- strsplit(sample, "R")[[1]][1]
      
      if(!grow_sat){#first option using t_vector as number of preyrations
        for(t in 1:t_vector[sample]){
          #x[, sample]<- x[,sample] +  t(mapply(function(x,eik)rbinom(n=1, size=x, prob=eik), x[, sample], qs$e_ik))
          x[, sample]<- x[,sample] +  t(mapply(function(x,eik)qbinom(p=runif(1), size=x, prob=eik, F), x[, sample], qs$e_ik))
          
        }
      }else{#second option using a while loop until we reach the number of cells 7.5e10
        while(sum(x[, sample])<=7.5e10){
          #x[, sample]<- x[,sample] +  t(mapply(function(x,eik)rbinom(n=1, size=x, prob=eik), x[, sample], qs$e_ik))
          x[, sample]<- x[,sample] +  t(mapply(function(x,eik)qbinom(p=runif(1), size=x, prob=eik, F), x[, sample], qs$e_ik))
          
        }
      }
      
      z[, sample]<- ceiling(t(mapply(function(m,phi)rnbinom(1, 1/phi, mu=m), x[, sample]*L_vector[sample], qs$phi_ikS)))
    }
    #z[is.na(z[,grep(paste0(bait,"R"), colnames(z))]),grep(paste0(bait,"R"), colnames(z))]<-0
    count_list[[bait]]<- z[, grep(paste0(bait,"R"), colnames(z))]
  }
  return(count_list)
}


simulate_prey_counts<-function(Q_N, Q_S, L_vector, t_vector, ns_gen=4, num_true_int, num_rep, num_baits, n_preys, 
                               true_int_thr, stickiness, M_0=3.84e9, low_conc=F, high_phi=F, overdispersed=F){
  library(tidyverse)
  #create database for sampling true int, sticky and non-interactors
  quantile_eik<- quantile(Q_S[Q_S$e_ik>0,"e_ik"], true_int_thr)
  Q_S_false_int<- Q_S[Q_S$e_ik< quantile_eik,]
  quantile_eik_sticky<- quantile(Q_S_false_int$e_ik, 0.99)
  Q_S_sticky<- Q_S_false_int[Q_S_false_int$e_ik> quantile_eik_sticky,]
  Q_S_sticky$interactor<- "S"
  Q_N_sim<- data.frame(prey= paste0("prey_", 1:n_preys),Q_N[sample(x=1:nrow(Q_N), size=n_preys, replace = T), 2:3])
  Q_N_sim$q_k<- Q_N_sim$q_k/sum(Q_N_sim$q_k)
  Q_S_sim<- simulate_non_interactors(Q_S, paste0('prey_', 1:n_preys), num_baits)
  Q_S_sim$interactor<- "N"
  
  sticky_preys<- paste0('prey_', sample(x=1:n_preys,size=round(stickiness*n_preys), replace = F))
  Q_S_sim[Q_S_sim$prey %in% sticky_preys, c("e_ik","phi_ikS","interactor")]<- Q_S_sticky[sample(x=1:nrow(Q_S_sticky),size=num_baits*round(stickiness*n_preys), replace = T), c("e_ik","phi_ikS","interactor")]
  
  bait= unlist(mapply(function(c,x)rep(c, each=x), c=names(num_true_int),x=num_true_int))
  
  Q_S_true<- Q_S[Q_S$e_ik>= quantile_eik,]
  true_int_index<- sample(x=1:nrow(Q_S_true),size=sum(num_true_int), replace = T)
  true_int<- data.frame(prey=sample(unique(Q_S_sim$prey)[-grep(paste(sticky_preys,collapse = "|"), unique(Q_S_sim$prey))], sum(num_true_int)), 
                        bait= unlist(mapply(function(c,x)rep(c, each=x), c=names(num_true_int),x=num_true_int)),
                        e_ik= Q_S_true[true_int_index,"e_ik"], phi_ikS= Q_S_true[true_int_index,"phi_ikS"],
                        interactor= rep("T", sum(num_true_int)), stringsAsFactors = F)
  for(bait in unique(true_int$bait)){
    Q_S_sim[Q_S_sim$prey %in% true_int[true_int$bait==bait, "prey"] & Q_S_sim$bait==bait, c("e_ik","phi_ikS","interactor")]<- true_int[true_int$bait==bait, c("e_ik","phi_ikS","interactor")]
  }
  
  count_list<- galton_walton_sim_table(Q_S_sim, Q_N_sim, t_vector, L_vector, num_preys=n_preys, num_rep, M_0, ns_gen, low_conc, true_int, high_phi, overdispersed)
  return(list(count_list, Q_S_sim, Q_N_sim))
}

#Now for fusion counts

simulate_fusion_counts<-function(count_list, Q_S_sim, Q_N_sim, num_rep, n_preys, u_kn, u_iks, p_kn, p_iks, num_baits, num_true_int){
  library(tidyverse)
  fusion_count_list <- list()
  u_kn<- sample(u_kn, size=n_preys,replace = T)
  p_kn<- sample(p_kn, size=n_preys,replace = T)
  Q_N_sim$u_kn<- u_kn
  Q_N_sim$p_kn<- p_kn
  qiks_true<- p_iks[p_iks>=quantile(p_iks, 0.95)]
  Q_S_sim$u_iks<-0
  Q_S_sim$p_iks<-0
  
  for (sample in names(count_list)){
    n_ns<- t(mapply(function(n)rbinom(n=num_rep, size=count_list[[sample]][,grep("N", colnames(count_list[[sample]]))][n,], prob=u_kn[n]), 1:nrow(count_list[[sample]])))
    np_ns<- t(mapply(function(n)rbinom(n=num_rep, size=n_ns[n,], prob=p_kn[n]), 1:nrow(n_ns)))
    
    u_iks_sample<- sample(u_iks, size=n_preys,replace = T)
    n_s<-t(mapply(function(n)rbinom(n=num_rep, size=count_list[[sample]][,-grep("N", colnames(count_list[[sample]]))][n,], prob=u_iks[n]), 1:nrow(count_list[[sample]])))
    
    p_iks_sample<- sample(p_iks[p_iks<=quantile(p_iks, 0.95)], size=n_preys,replace = T)
    p_iks_sample[as.numeric(sapply(Q_S_sim[Q_S_sim$bait==sample & Q_S_sim$interactor=="T", "prey"], FUN=function(x)strsplit(x, "_")[[1]][2]))]<- sample(qiks_true, size=num_true_int[[sample]], replace = T)
    np_s<- t(mapply(function(n)rbinom(n=num_rep, size=n_s[n,], prob=p_iks_sample[n]), 1:nrow(n_ns)))
    Q_S_sim[Q_S_sim$bait==sample, "u_iks"]<-u_iks_sample
    Q_S_sim[Q_S_sim$bait==sample, "p_iks"]<-p_iks_sample
    # put everything together reps in rows 
    s_preys<- Q_S_sim[Q_S_sim$bait==sample, "prey"]
    
    fusion_table_np_s<- data.frame(Transcript_id= paste0(s_preys, ".1"),np_s, stringsAsFactors = F)
    colnames(fusion_table_np_s)<- c("Transcript_id", paste0(rep(sample, each=num_rep),"R",1:num_rep, "S"))
    fusion_table_np_s<- gather(fusion_table_np_s, Replicate, num_fusion_reads_in_frame_selected_sample, 2:(num_rep+1))
    
    fusion_table_n_s<- data.frame(Transcript_id= paste0(s_preys, ".1"),n_s, stringsAsFactors = F)
    colnames(fusion_table_n_s)<- c("Transcript_id", paste0(rep(sample, each=num_rep),"R",1:num_rep, "S"))
    fusion_table_n_s<- gather(fusion_table_n_s, Replicate, num_fusion_reads_selected_sample, 2:(num_rep+1))
    
    fusion_table_np_ns<- data.frame(Transcript_id= paste0(s_preys, ".1"),np_ns, stringsAsFactors = F)
    colnames(fusion_table_np_ns)<- c("Transcript_id", paste0(rep(sample, each=num_rep),"R",1:num_rep, "S"))
    fusion_table_np_ns<- gather(fusion_table_np_ns, Replicate, num_fusion_reads_in_frame_background_sample, 2:(num_rep+1))
    
    fusion_table_n_ns<- data.frame(Transcript_id= paste0(s_preys, ".1"),n_ns, stringsAsFactors = F)
    colnames(fusion_table_n_ns)<- c("Transcript_id", paste0(rep(sample, each=num_rep),"R",1:num_rep, "S"))
    fusion_table_n_ns<- gather(fusion_table_n_ns, Replicate, num_fusion_reads_background_sample, 2:(num_rep+1))
    
    ns_table<- cbind(fusion_table_np_ns, fusion_table_n_ns[,3])
    s_table<- cbind(fusion_table_np_s, fusion_table_n_s[,3])
    
    fusion_count_list[[sample]]<- merge(s_table, ns_table, by=c("Transcript_id","Replicate"), all=T)
    fusion_count_list[[sample]][is.na(fusion_count_list[[sample]])] <- 0
    colnames(fusion_count_list[[sample]])<- c("Transcript_id","Replicate", "num_fusion_reads_in_frame_selected_sample",
                                              "num_fusion_reads_selected_sample", "num_fusion_reads_in_frame_background_sample",
                                              "num_fusion_reads_background_sample")
    
  }
  return(list(fusion_count_list, Q_S_sim, Q_N_sim))
}


# lets save the files in the required format
save_files<- function(QS, QN, total_count_list, fusion_count_list, output_directory){
  dir.create(output_directory)
  for (sample in names(total_count_list)) {
    write.table(total_count_list[[sample]], paste0(output_directory,sample,"_salmon_counts.matrix"), row.names = T, col.names = F)
    write.csv(fusion_count_list[[sample]], paste0(output_directory,sample,"_final_report.csv"))
  }
  saveRDS(total_count_list, paste0(output_directory, "total_count_list.RDS"))
  saveRDS(fusion_count_list, paste0(output_directory, "fusion_count_list.RDS"))
  saveRDS(QS, paste0(output_directory, "QS.RDS"))
  saveRDS(QN, paste0(output_directory, "QN.RDS"))
}


#run the scores to see results
# check the performance of the scores
performance<- function(input_files_location, out_dir){
  list_true_int_preys <- readRDS(paste0(input_files_location, "QS.RDS"))[,c("bait","prey","interactor")]
  list_true_int_preys$interaction<- list_true_int_preys$interactor=="T"
  #colnames(list_true_int_preys)<- c("bait","prey","interaction")
  list_true_int_preys<- list_true_int_preys[,c("bait","prey","interaction")]
  total_scores <- readRDS(paste0(input_files_location,"output/total_scores.RDS"))
  total_scores<- merge(total_scores, list_true_int_preys, by=c("bait","prey"), all.x=T)
  pdf(paste0(out_dir,"performance_plots.pdf"), width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
  library(plotROC)
  library(tidyverse)
  basicplot1<- ggplot(total_scores, aes(m = total_bait_enrich_score, d = interaction))+ geom_roc(n.cuts=2,labels=T, labelsize = 20) 
  print(basicplot1 + style_roc(theme = theme_grey)+ geom_rocci(fill="violet", labelsize = 20) + theme_grey(base_size = 22)+
          annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot1)$AUC, 2)), size = 20))
  
  basicplot2<- ggplot(total_scores, aes(m = total_spec_score, d = interaction))+ geom_roc(n.cuts=2,labels=T, labelsize = 20) 
  print(basicplot2 + style_roc(theme = theme_grey) + geom_rocci(fill="cyan", labelsize = 20) + theme_grey(base_size = 22)+
          annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot2)$AUC, 2)), size = 20))
  
  basicplot3<- ggplot(total_scores, aes(m = total_max_freq_score, d = interaction))+ geom_roc(n.cuts=2,labels=T, labelsize = 20)
  print(basicplot3 + style_roc(theme = theme_grey) + geom_rocci(fill="green", labelsize = 20) + theme_grey(base_size = 22)+
          annotate("text", x = .75, y = .25, label = paste("AUC =", round(calc_auc(basicplot3)$AUC, 2)), size = 20))
  
  
  #dev.off()
  #PR curve
  library(PRROC)
  #pdf(paste0(out_dir,"PR_curves_newScores_nf.pdf"), width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
  
  pr_2<-pr.curve(scores.class0 = total_scores[!is.na(total_scores$interaction) & !is.na(total_scores$total_bait_enrich_score), "total_bait_enrich_score"], 
                 weights.class0 = total_scores[!is.na(total_scores$interaction)& !is.na(total_scores$total_bait_enrich_score), "interaction"], curve=T)
  plot(pr_2)
  
  pr_3<-pr.curve(scores.class0 = total_scores[!is.na(total_scores$interaction)& !is.na(total_scores$total_spec_score), "total_spec_score"], 
                 weights.class0 = total_scores[!is.na(total_scores$interaction)& !is.na(total_scores$total_spec_score), "interaction"], curve=T)
  plot(pr_3)
  
  pr_1<-pr.curve(scores.class0 = total_scores[!is.na(total_scores$interaction)&!is.na(total_scores$total_max_freq_score), "total_max_freq_score"], 
                 weights.class0 = total_scores[!is.na(total_scores$interaction)&!is.na(total_scores$total_max_freq_score), "interaction"], curve=T)
  plot(pr_1)
  
  #dev.off()
  #put auc in the scenario table
  auc<- data.frame(ROC_AUC_Enrichment=round(calc_auc(basicplot1)$AUC, 2), ROC_AUC_Specificity=round(calc_auc(basicplot2)$AUC, 2), ROC_AUC_in_frame=round(calc_auc(basicplot3)$AUC, 2),
                   PR_AUC_Enrichment=pr_2$auc.integral, PR_AUC_Specificity=pr_3$auc.integral, PR_AUC_in_frame=pr_1$auc.integral, stringsAsFactors = F)
  
  #other plot
  #library(tidyverse)
  #pdf(paste0(out_dir,"density_scores_newScores.pdf"), width = 15, height = 15, fonts = "ArialMT", pointsize = 14)
  library(tidyverse)
  print(ggplot(total_scores[!is.na(total_scores$interaction), ], aes(total_spec_score, fill=as.factor(interaction))) + 
          geom_density(alpha=0.5) +scale_fill_manual(values = c("darkgreen", "violet"))+ theme_grey(base_size = 22))
  
  print(ggplot(total_scores[!is.na(total_scores$interaction), ], aes(total_bait_enrich_score, fill=as.factor(interaction))) + 
          geom_density(alpha=0.5) +scale_fill_manual(values = c("red", "blue"))+ theme_grey(base_size = 22))
  
  print(ggplot(total_scores[!is.na(total_scores$interaction), ], aes(total_max_freq_score, fill=as.factor(interaction))) + 
          geom_density(alpha=0.5) +scale_fill_manual(values = c("darkgreen", "orange"))+ theme_grey(base_size = 22))
  
  dev.off()
  return(auc)
}

#preload
Q_N<- readRDS("Q_N.RDS")
Q_N[Q_N<0]=0
Q_N<-Q_N[1:25575,2:4]
Q_S<- readRDS("Q_S.RDS")
Q_S[Q_S<0]=0
L_vector<- readRDS("L_vector.RDS")
t_vector_s<- readRDS("ts_vector.RDS")[names(L_vector)[-grep("N", names(L_vector))]]
t_vector_n<-rep(4,30)
names(t_vector_n)<-gsub('.{1}$', 'N', names(t_vector_s))
U<- readRDS("U.RDS")
P<- readRDS("P.RDS")
u_kn<- U[U$key=="u_kn","value"]
u_iks<- U[U$key=="u_iks","value"]
p_kn<- P[P$key=="pi_kn","value"]
p_iks<- P[P$key=="pi_iks","value"]


#I created a table with all the parameters to run the scores and save the AUC values
scenarios<- read.csv("simulation_scenarios.csv", stringsAsFactors = F)
dir.create("scenarios/")

for (i in 1:nrow(scenarios)) { #nrow(scenarios)
  num_baits=scenarios$Baits[i]
  sample_names<- paste0("bait", 1:num_baits)
  num_rep<-scenarios$Replicates[i]
  n_preys<- scenarios$Preys[i]
  size_true<- ceiling(runif(length(sample_names), n_preys*0.0004, n_preys*0.001))
  names(size_true)<- sample_names
  t_vector_1<- sample(t_vector_s, size= num_rep*num_baits, replace = T)
  names(t_vector_1)<- paste0(rep(sample_names, each=num_rep),"R",1:num_rep, "S")
  L_vector_1<- sample(L_vector, size= 2*num_rep*num_baits, replace = T)
  names(L_vector_1)<- c(paste0(rep(sample_names, each=num_rep),"R",1:num_rep, "S"),
                        paste0(rep(sample_names, each=num_rep),"R",1:num_rep, "N"))
  stickiness<- scenarios$Stickiness[i]
  true_int_thr<-scenarios$TI_strength[i]
  low_conc<- scenarios$TI_t0[i]=="Low"
  #high_phi<- scenarios$Overdispersion[i]=="High"
  overdispersed<-scenarios$Overdispersion[i]=="High"
  simulation_1<- simulate_prey_counts(Q_N, Q_S, L_vector_1, t_vector_1, ns_gen=4, size_true, num_rep, num_baits, n_preys, true_int_thr, stickiness, M_0=3.84e9, low_conc, high_phi=F, overdispersed)
  #simulation_1<- simulate_prey_counts(Q_N, Q_S, L_vector_1, t_vector_1, ns_gen=4, size_true, num_rep, num_baits, n_preys, true_int_thr, stickiness, M_0=3.84e9, low_conc, high_phi=F, overdispersed = T)
  
  f_simulation_1_s<- simulate_fusion_counts(simulation_1[[1]], simulation_1[[2]], simulation_1[[3]], num_rep, n_preys, u_kn, u_iks, p_kn, p_iks, num_baits, size_true)
                                            
  #input_files_location<-paste0("scenarios/", scenarios$Scenario[i],"/")
  input_files_location<-paste0("scenarios/", scenarios$Scenario[i],"/")
  
  total_count_list<- simulation_1[[1]]
  fusion_count_list<- f_simulation_1_s[[1]]
  QS<- f_simulation_1_s[[2]]
  QN<- f_simulation_1_s[[3]]
  save_files(QS, QN, total_count_list, fusion_count_list, input_files_location)
  
  #run scores
  location_folders<- "scores/"
  out_dir<- paste0(input_files_location, "output/")
  p_val_thr<- 0.5
  p_val_thr_s<-0.01
  fc_thr<-2
  num_replicates<- rep(num_rep,length(sample_names))
  norm_plots<-T
  source(paste0(location_folders,"run_scores_srcFun.R"))
  run_scores(sample_names, num_replicates, input_files_location, location_folders, out_dir, p_val_thr, p_val_thr_s, fc_thr, norm_plots)
  #see performance
  auc<- performance(input_files_location, out_dir)
  scenarios[i, 9:14]<- auc
  
}

write.csv(scenarios, "simulation_scenarios_computed.csv")



