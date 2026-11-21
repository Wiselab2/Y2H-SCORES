#get interacting sequence script
#  interacting_fragments<- get_interaction_files(fasta_list, sample_names, final_reports_list, output_location)

get_interaction_files<- function(fasta_list, total_scores, sample_names, final_reports_list, output_location, thr=0){
  require(tidyverse)
  require(Biostrings)
  
  getmode <- function(v) {
    uniqv <- sort(unique(v), decreasing = T)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  fastas<-list()
  for(num in 1:length(sample_names)){
    sample <- sample_names[num]
    final_report<- read.csv(final_reports_list[num],header = TRUE)
    FL<- readDNAStringSet(fasta_list[num])
    total_scores_sample<-total_scores[total_scores$bait==sample,]
    
    final_report$junction_read_start_mode<-NA
    final_report[final_report$starting_location_of_in_frame_junction_reads!="", "junction_read_start_mode"]<-sapply(sapply(final_report[final_report$starting_location_of_in_frame_junction_reads!="","starting_location_of_in_frame_junction_reads"], FUN=function(x) strsplit(x,";")),
                                                                                                                    FUN=function(x) getmode(as.numeric(unlist(sapply(x, FUN=function(y) strsplit(y,":")[[1]][2])))))
    final_report$max_stop_5UTR<--Inf
    final_report[final_report$location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS!="","max_stop_5UTR"]<-sapply(final_report[final_report$location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS!="","location_of_STOP_codons_in_5_prime_UTR_region_in_frame_with_CDS"], FUN=function(x) max(as.numeric(unlist(strsplit(x,";")))))
    
    
    final_report$start_interaction<- ifelse(!is.na(final_report$junction_read_start_mode) & final_report$junction_read_start_mode >final_report$max_stop_5UTR, 
                                            final_report$junction_read_start_mode, final_report$CDS_start)
    
    
    final_report$extraction_coords<- paste0(paste(final_report$Transcript_id,final_report$Replicate, sep = "_"),":", 
                                            final_report$start_interaction, "-",final_report$CDS_end)
    
    a<- sapply(final_report[final_report$starting_location_of_in_frame_junction_reads!="","starting_location_of_in_frame_junction_reads"], FUN=function(x) strsplit(x,";"))
    names(a)<- paste0(final_report[final_report$starting_location_of_in_frame_junction_reads!="","Transcript_id"], "_", final_report[final_report$starting_location_of_in_frame_junction_reads!="","Replicate"])
    b<- sapply(a, FUN=function(x) sapply(x, FUN=function(y) as.numeric(strsplit(y,":")[[1]][2])))
    n<- sapply(a, FUN=function(x) sapply(x, FUN=function(y) strsplit(y,":")[[1]][1]))
    
    final_report$add_seq<- ""
    rows<-as.numeric(row.names(final_report[final_report$starting_location_of_in_frame_junction_reads!="",]))
    for(i in 1:length(a)){
      if(final_report$start_interaction[rows[i]]==final_report$junction_read_start_mode[rows[i]]){
        name<- names(a)[i]
        hit<-n[[name]][b[[name]]==final_report$junction_read_start_mode[rows[i]]]
        #if(length(grep("N", hit))>0 & length(hit)>1){
        if(length(grep("N", hit)) < length(hit) & length(grep("N", hit))>0){ 
          c<-table(hit[-grep("N", hit)])
        }else{
          c<-table(hit)
        }
        final_report$add_seq[rows[i]]<- names(c[which(c==max(c))])[1]
      }
    }
    
    #final_report_filtered<- final_report[final_report$Gene %in% total_scores_sample[total_scores_sample$Sum_scores>thr, "prey"],]
    
    final_report_filtered<- final_report[final_report$Gene %in% total_scores_sample[total_scores_sample$Sum_scores>thr, "prey"] &
                                          final_report$num_fusion_reads_selected_sample>0 & final_report$num_fusion_reads_background_sample>0,]
    
    
    fasta_name<- paste0(final_report_filtered$Transcript_id, "_", final_report_filtered$Replicate)
    final_report_filtered$fasta_name<- paste0(final_report_filtered$Transcript_id, "_", final_report_filtered$Replicate)
    #names(FL)
    FL_filtered<- FL[fasta_name[fasta_name %in% names(FL)]]
    add_seq<- final_report_filtered$add_seq
    names(add_seq)<- paste0(final_report_filtered$Transcript_id, "_", final_report_filtered$Replicate)
    add_seq_DNA <- DNAStringSet(add_seq)
    
    interaction_frags<- FL_filtered
    interaction_report_filtered<- final_report_filtered[final_report_filtered$fasta_name %in% names(interaction_frags), ]
    for(i in 1:nrow(interaction_report_filtered)){
      name<-paste0(interaction_report_filtered$Transcript_id[i], "_", interaction_report_filtered$Replicate[i])
      interaction_frags[name]<- xscat(add_seq[name], subseq(FL_filtered[name], start=interaction_report_filtered$start_interaction[i], end=interaction_report_filtered$CDS_end[i]))
    }
    fastas[[sample]]<- interaction_frags
    #write.csv(final_report, paste0("~/iowa_state/lab/Y2H-Seq effectors/results/SE_", effector,"/SE_", effector,"_final_report_position_fusion_seq.csv"), row.names = FALSE)
    #saveRDS(in_frame_results_sum, paste0(output_location, "in_frame_results_sum.RDS"))
    write.csv(final_report_filtered, paste0(output_location, sample,"_final_report_position_fusion_seq.csv"), row.names = FALSE)
    writeXStringSet(interaction_frags, paste0(output_location, sample,"_interaction_fragments_inframe.fasta"))
    writeXStringSet(FL_filtered, paste0(output_location, sample,"_transcriptome_file_full_length.fasta"))
  }
  return(fastas)
}

