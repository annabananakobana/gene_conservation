#author=Alana Alexander
# load necessarily libraries
# install.packages("tidyverse")
library(tidyverse)
# library(devtools)
# devtools::install_github("mhahsler/rBLAST")
library(rBLAST)

# Adding the location of the blast exectuables
# to R's path
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/nesi/CS400_centos7_bdw/BLAST/2.10.0-GCC-9.2.0/bin", sep= .Platform$path.sep))

# Getting a list of genes, and excluding no_length_filter
list_files <- list.files("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files")

list_files <- list_files[which(list_files!="discard")]

sequence_location_record <- NULL

# Stepping through each gene
for (i in list_files) {
  print(paste("We are up to gene ",i,sep=""))
  # Getting a list of the files we are interested in
  fasta_files <- list.files(paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,sep=""), pattern="*_labelled.fasta")
  maf_file <- list.files(paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,"/results_5_20_75",sep=""),pattern="*m80.maf")
  
  # Reading in the maf file
  read_in_maf_file <- readLines(paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,"/results_5_20_75/",maf_file,sep=""))
  
  # Numbering matches in the same area of the genome sequentially
  row_numbers <- seq(1,length(read_in_maf_file))[grep("^s ",read_in_maf_file)]
  match_numbers <- rep(NA,length(row_numbers))
  match_numbers[1] <- 1
  for (j in 2:length(match_numbers)) {
    if (row_numbers[j]==(row_numbers[j-1] + 1 )) {
      match_numbers[j] <- match_numbers[j-1]
    } else {
      match_numbers[j] <- match_numbers[j-1] + 1
    }
  }
  
  # Taking just the lines with matches on them
  read_in_maf_file <- read_in_maf_file[grep("^s ",read_in_maf_file)]
  
  # Converting into a table
  number_of_matches <- length(read_in_maf_file)
  read_in_maf_file <- matrix(unlist(strsplit(read_in_maf_file,"\\s+")),nrow=number_of_matches,byrow = TRUE)
  # Keeping just the columns we are interested in
  read_in_maf_file <- read_in_maf_file[,c(2,7)]
  
  # Creating a column where the gaps have been stripped out for blasting
  # And adding in our match numbers column
  # Adding in a column for beginning and end
  # of matching
  read_in_maf_file <- cbind(read_in_maf_file,gsub("-","",read_in_maf_file[,2]),match_numbers,NA,NA)
  
  # Making a blast database for each of the fasta files
  for (j in fasta_files) {
    makeblastdb(paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,"/",j,sep=""), dbtype = "nucl", args="")
  }
  
  # Because maf file names do not match up to sequence name
  # constructing a key based on the sequence headers of
  # each of the fasta files
  sequence_header_key <- matrix(NA,ncol=3,nrow=length(fasta_files))
  for (j in 1:length(fasta_files)) {
    sequence_header_key[j,1] <- fasta_files[j]
    sequence_header_key[j,2] <- gsub(":.*","",gsub(">","",readLines(paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,"/",fasta_files[j],sep=""))[1]))
  }
  
  sequence_header_key <- sequence_header_key[order(sequence_header_key[,1]),]
  
  sequence_header_key <- cbind(sequence_header_key,c("red","orange","yellow3","yellowgreen","green3","green4","deepskyblue","dodgerblue4","darkviolet"))
  
  # Removing any lines for genomes where no sequence is present
  # in the fasta files
  if(any(is.na(sequence_header_key[,2]))) {
    sequence_header_key <- sequence_header_key[-which(is.na(sequence_header_key[,2])),]
  }
  
  # Additional rows to record
  additional_rows <- NULL
  
  # break variable if mismatch between maf and fasta
  break_on_gene <- NULL
  
  # Running blast to get beginning and end of each match
  for (j in 1:dim(read_in_maf_file)[1]) {
    
    # Finding out what fasta file we should blast against
    k <- which(sequence_header_key[,2]==read_in_maf_file[j,1])
    
    if(length(k)==0) {
      print(paste(read_in_maf_file[j,1]," sequence in the maf file",sep=""))
      print("does not exist among the fasta files. Perhaps fasta files are")
      print("from a different version of reference? Please ensure fasta")
      print("file version match maf alignments and try again")
      break_on_gene <- i
      break
    }
    
    sequence_header_key[k,3] <- read_in_maf_file[j,1]
    
    # Loading in database based on this name
    bl <- blast(db = paste("/nesi/nobackup/uoo02820/gene_conservation_multiz/assembly_databases/genomes/final_fasta_files/",i,"/",sequence_header_key[k,1],sep=""))
    
    # blasting to the database
    cl <- predict(bl,DNAStringSet(read_in_maf_file[j,3]))
    
    # Only taking the perfect matchees 
    if (dim(cl)[1] > 1) {
      if(any(cl$Mismatches>0)) {
        if(length(which(cl$Mismatches>0))<dim(cl)[1]) {
          cl <- cl[-which(cl$Mismatches>0),]
        }
      }  
      if(length(which(cl$Alignment.Length==min(cl$Alignment.Length)))<dim(cl)[1]) {
        cl <- cl[-which(cl$Alignment.Length==min(cl$Alignment.Length)),]  
      }
    }
    
    # And if there is still more than one row...
    if (dim(cl)[1] > 1) {
      # Getting the start coordinates
      read_in_maf_file[j,5] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[1,2]))) + as.numeric(cl[1,9]) - 1
      
      # Getting the end coordinates
      read_in_maf_file[j,6] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[1,2]))) + as.numeric(cl[1,10]) - 1
      
      for (k in 2:dim(cl)[1]) {
        temp_additional_rows <- read_in_maf_file[j,]
        # Getting the start coordinates
        temp_additional_rows[5] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[k,2]))) + as.numeric(cl[k,9]) - 1
        
        # Getting the end coordinates
        temp_additional_rows[6] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[k,2]))) + as.numeric(cl[k,10]) - 1
        
        # These additional rows being recorded
        additional_rows <- rbind(additional_rows,temp_additional_rows)
      }
    } else {  
      # Getting the start coordinates
      read_in_maf_file[j,5] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[2]))) + as.numeric(cl[9]) - 1
      
      # Getting the end coordinates
      read_in_maf_file[j,6] <- as.numeric(gsub("-.*","",gsub(".*:","",cl[2]))) + as.numeric(cl[10]) - 1
    }
    
  }
  
  if (!is.null(break_on_gene)) {
    next
  }
  
  read_in_maf_file <- rbind(read_in_maf_file,additional_rows)
  
  # Converting to a tibble
  read_in_maf_file <- as_tibble(read_in_maf_file)
  read_in_maf_file <- read_in_maf_file  %>% mutate_at(vars(match_numbers:V6),funs(as.numeric))
  names(read_in_maf_file) <- c("genome","seq_w_gaps","seq_wo_gaps","match_numbers","start","end")
  
  # Getting the number of matching stretches
  num_matches <- length(unique(read_in_maf_file$match_numbers))
  
  # Identifying the benchmark genome i.e. the genome
  # that is present in all interactions
  # breaking from the loop and identifying what loci
  # if no genome present in all comparisons
  benchmark <- read_in_maf_file %>% group_by(genome) %>% count(match_numbers) %>% count() %>% 
    filter(n==num_matches) %>% select(genome) %>% as.matrix()
  
  benchmark <- benchmark[1,1]
  
  if (length(benchmark)==0) {
    print(paste(i, " does not have a single genome present in all matches. Examine manually",sep=""))
    break
  }
  
  desired_order <- ((read_in_maf_file %>% mutate(totalmatch=abs(end-start)) %>% group_by(genome) %>% summarize(summatches=sum(totalmatch)) %>% arrange(desc(summatches)) %>% select(genome) %>% as.matrix())[,1])
  
  max_height <- length(desired_order)
  
  # Getting the data together to plot base layers
  min_x <- min(read_in_maf_file %>% filter(genome==benchmark) %>% select(start,end))
  max_x <- max(read_in_maf_file %>% filter(genome==benchmark) %>% select(start,end))
  
  line_plot_data <- matrix(NA,ncol=5,nrow=length(desired_order))
  line_plot_data[,1] <- desired_order
  line_plot_data[,2] <- min_x
  line_plot_data[,3] <- max_x
  line_plot_data[,4] <- rev(seq(2,length(desired_order)*3-1,3))
  line_plot_data[,5] <- sequence_header_key[match(desired_order,sequence_header_key[,3]),4]
  
  line_plot_data <- as_tibble(line_plot_data)
  
  names(line_plot_data) <- c("genome","min_x","max_x","ys","colours")
  
  line_plot_data <- line_plot_data %>% mutate_at(vars(min_x:ys),funs(as.numeric))
  
  line_plot_data$genome <- factor(line_plot_data$genome,levels=desired_order)
  
  # Plotting base layers first
  output_plot <- ggplot() + 
    geom_segment(data=line_plot_data,aes(x=min_x,xend=max_x,y=ys,yend=ys,color=genome),size=3) +
    scale_colour_manual(values=as.matrix(line_plot_data[,5])[,1])
  
  # Getting data together to plot conserved chunks
  match_key <- read_in_maf_file %>% filter(genome==benchmark) %>% select(genome,match_numbers,start,end)
  names(match_key) <- c("genome","match_numbers","ref_start","ref_end")
  
  conserved_plot_data <- read_in_maf_file %>% select(genome,match_numbers,start,end)
  
  conserved_plot_data <- full_join(conserved_plot_data,match_key,"match_numbers")
  
  conserved_plot_data <- conserved_plot_data %>% select(genome.x,match_numbers,ref_start,ref_end) %>% distinct()
  
  conserved_plot_data <- conserved_plot_data %>% rowwise() %>%  mutate(smallest=min(ref_start,ref_end),biggest=max(ref_start,ref_end)) %>% select(genome.x,match_numbers,smallest,biggest)
  
  tempsequence_location_record <- full_join(conserved_plot_data,read_in_maf_file,"match_numbers") %>% select(genome,match_numbers,smallest,biggest,seq_wo_gaps,start,end) %>% distinct()
  tempsequence_location_record <- cbind(i,tempsequence_location_record)
  sequence_location_record <- rbind(sequence_location_record,tempsequence_location_record)
  
  conserved_plot_data <- cbind(conserved_plot_data,NA,NA)
  
  names(conserved_plot_data) <- c("genome","match_numbers","min_x","max_x","min_y","max_y")
  
  conserved_plot_data$min_y <- line_plot_data$ys[match(conserved_plot_data$genome,line_plot_data$genome)]-1
  
  conserved_plot_data$max_y <- line_plot_data$ys[match(conserved_plot_data$genome,line_plot_data$genome)]+1
  
  conserved_plot_data$genome <- factor(conserved_plot_data$genome,levels=desired_order)  
  
  output_plot + geom_rect(data=conserved_plot_data,inherit.aes=FALSE,aes(xmin=min_x,xmax=max_x,ymin=min_y,ymax=max_y,group=match_numbers,fill=genome),color="black") + 
    scale_fill_manual(values=as.matrix(line_plot_data[,5])[,1]) +
    theme_bw(base_size=18) +
    theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
    theme(axis.title=element_text(size=20,face="bold")) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
    scale_x_continuous(name=paste("Distance (bp) relative to ",benchmark,sep="")) +
    labs(fill="Genome") +
    guides(color=FALSE)
  
  ggsave(filename=paste(i,"_conserved_blocks.pdf",sep=""),plot = last_plot(),width=297,height=210,units="mm")
  
  
}

names(sequence_location_record) <- c("Gene","genome","match_numbers","start_relative_to_benchmark","end_relative_to_benchmark","seq_wo_gaps","start","end")

write_delim(sequence_location_record,"match_sequence_summary.txt")
