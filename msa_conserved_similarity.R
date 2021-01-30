#author=Alana Alexander
# load necessarily libraries
# install.packages("tidyverse")
library(tidyverse)
# library(devtools)
# devtools::install_github("mhahsler/rBLAST")
library(rBLAST)
# devtools::install_github("UBod/msa")
library(msa)

# Adding the location of the blast exectuables
# to R's path
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/nesi/CS400_centos7_bdw/BLAST/2.10.0-GCC-9.2.0/bin", sep= .Platform$path.sep))

# Reading in our data created in the visualise_conservation loop
sequence_location_record <- read_delim("match_sequence_summary.txt",delim = " ")

output <- NULL

# Stepping through each gene
for (i in unique(sequence_location_record$Gene)) {
  print(paste("We are up to gene ",i,sep=""))
  
  gene_sequence_location_record <- sequence_location_record %>% filter(Gene==i) %>% arrange()
  
  # "Grouping" matches that overlap in genomic coordinates
  unique_coords <- gene_sequence_location_record %>% select(start_relative_to_benchmark,end_relative_to_benchmark) %>% distinct() %>% arrange(start_relative_to_benchmark)
  unique_coords <- cbind(unique_coords,NA)
  names(unique_coords)[3] <- "combined_match"
  
  unique_coords[1,3] <- 1
  
  if (dim(unique_coords)[1]>1) {
    for (j in 2:dim(unique_coords)[1]) {
      if (unique_coords[j,1] < unique_coords[(j-1),2]) {
        unique_coords[j,3] <- unique_coords[(j-1),3]
      } else {
        unique_coords[j,3] <- unique_coords[(j-1),3] + 1
      }
    }
  }
  
  gene_sequence_location_record <- inner_join(gene_sequence_location_record,unique_coords,c("start_relative_to_benchmark","end_relative_to_benchmark"))
  
  # Creating a consensus sequence for each species represented multiple times in each combined_match
  # For each genome
  for (j in unique(gene_sequence_location_record$genome)) {
    tempspecies <- gene_sequence_location_record %>% filter(genome==j) %>% arrange(combined_match,match_numbers)
    
    # Ensuring all matches are given 5' to 3'
    if (tempspecies$start[1] < tempspecies$end[1]) {
      tempspecies <- tempspecies %>% mutate(start_coords=start,end_coords=end)
    } else {
      tempspecies <- tempspecies %>% mutate(start_coords=end,end_coords=start)
    }
    
    # For each of the combined matches for the species we are working on
    for (k in unique(tempspecies$combined_match)) {
      tempspecies_match <- tempspecies %>% filter(combined_match==k) %>% arrange(start_coords)
      
      # If there multiple sequences for this species and combined match
      if (!is.null(dim(tempspecies_match)[1]) & dim(tempspecies_match)[1]>1) {
        # Stepping through the sequences and combining them
        for (m in 2:dim(tempspecies_match)[1]) {
          # If the ranges of the sequence overlap...
          if (tempspecies_match$start_coords[m]>=tempspecies_match$start_coords[(m-1)] & tempspecies_match$start_coords[m]<tempspecies_match$end_coords[(m-1)]) {
              # Creating a consensus sequence
              tempsequence <- msaConsensusSequence(msa(DNAStringSet(tempspecies_match$seq_wo_gaps[c((m-1),m)])),type="upperlower",thresh=c(0,0))
              tempsequence <- toupper(gsub("?","N",gsub("-","",tempsequence),fixed=TRUE))
          } else {
              # If they don't overlap, concatenating strings (because there is a gap of unknown size between them)
              if (tempspecies_match$start_relative_to_benchmark[(m-1)] <= tempspecies_match$start_relative_to_benchmark[m]) {
                tempsequence <- toupper(paste(tempspecies_match$seq_wo_gaps[c(m,(m-1))],collapse=""))
              } else {
                tempsequence <- toupper(paste(tempspecies_match$seq_wo_gaps[c((m-1),m)],collapse=""))
              }
          }
          tempspecies_match$seq_wo_gaps[1:m] <- tempsequence          
        }
        # Replace the consensus sequence for this species for this match, wherever it occurs
        gene_sequence_location_record$seq_wo_gaps[which(gene_sequence_location_record$genome==j & gene_sequence_location_record$combined_match==k)] <- tempsequence 
      }
    }
  }
  
  # Obtaining one row per species per combined_match
  gene_sequence_location_record <- gene_sequence_location_record %>% select(-match_numbers,-start_relative_to_benchmark,-end_relative_to_benchmark,-start,-end) %>% distinct() %>% arrange(combined_match)
  
  # For each of the combined_matches
  for (j in unique(gene_sequence_location_record$combined_match)) {
    # obtaining the sequences
    tempsequence <- gene_sequence_location_record$seq_wo_gaps[which(gene_sequence_location_record$combined_match==j)]
    
    seqence_alignment <- msa(DNAStringSet(tempsequence))
    
    alignment_count <- as_tibble(t(consensusMatrix(seqence_alignment)))

    prop_ident <- nrow(alignment_count %>% filter_all(any_vars(. == length(tempsequence))))/nrow(alignment_count)
    
    comparisons <- 0
    pairwise_similarity <- 0
      
    for (k in 1:(length(tempsequence)-1)) {
      for (m in (k+1):length(tempsequence)) {
        comparisons <- comparisons + 1
        
        pairwisesequences <- tempsequence[c(k,m)]
        pairwisesequence_alignment <- msa(DNAStringSet(pairwisesequences))
        
        tempsequences_count <- as_tibble(t(consensusMatrix(pairwisesequence_alignment)))
        
        pairwise_similarity <- pairwise_similarity + nrow(tempsequences_count %>% filter_all(any_vars(. == 2)))/nrow(tempsequences_count)
      }
    }
    
    output <- rbind(output,c(i,unique(gene_sequence_location_record$combined_match)[j],dim(alignment_count)[1],length(tempsequence),
                             prop_ident,(pairwise_similarity/comparisons)))
  }
  
}  

output <- as_tibble(output)
    
names(output) <- c("Gene","Conserved_block","Length_bp","No_of_indivs","Prop_identical","Av_pairwise_sim")
    
output <- output %>% mutate_at(vars(Conserved_block:Av_pairwise_sim),funs(as.numeric))
  
# ggplot() + geom_point(output,mapping=aes(x=No_of_indivs,y=Av_pairwise_sim,color=Gene))

write_delim(output,"match_similarity_for_conserved_blocks.txt")
 
