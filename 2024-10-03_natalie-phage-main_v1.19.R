library(gsheet)      
library(kableExtra)  
library(staplr)      
# library(randomForest) # For random forest
# library(rpart) # For decision tree (recursive partitioning)

################################################################################
# Define Parameters
################################################################################

minAmpliconLength <- 5    # Minimum sequence length 
maxAmpliconLength <- 1000000000    # Maximum sequence length
minCoverage <- 5           # Minimum coverage for a variant to be called
minAF <- 2/3                # Minimum allele frequency for a variant to be called
minMQ <- 30                 # Minimum mapping quality for a variant to be called 
minQUAL <- 0                # Minimum base quality for a variant to be called
QCThreshold <- 95           # Minimum percent 50x coverage for a specimen to have a consensus sequence made.
file_met <- '2024-09-19_natalie-phage-metadata_v3.csv' # Associated metadata
phage_min_len <- 50 # Minimum length for it to be called a phage region of interest

file_ref_out <- '2024-09-19_reference_combined.csv'


################################################################################
# Define Functions
################################################################################

# Generates a consensus sequence.
make_consensus <- function(i, minAF, minCoverage, ref){
  vsp <- ref$sample_id[i]
  
  # Determine if there are enough reads to make a consensus.
  if(file.size(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv')) == 0){
    ref$error[i] <- paste0(ref$error[i],'|no consensus made - too few reads|')
    print('|no consensus made - too few reads|')
  }
  if(file.size(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv')) > 0){
    
    # Open the samtools pileup
    pil1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_pileup.tsv'), sep='\t',header=FALSE)
    colnames(pil1) <- c('CHROM','POS','REF','COVERAGE','READ','QAUL')
    
    # Determine the occurrence of base per position.
    pil2 <- pil1
    pil2$Consensus <- NA
    pil2$A_count <- 0
    pil2$C_count <- 0
    pil2$G_count <- 0
    pil2$T_count <- 0
    # Reorder columns
    pil2 <- pil2[ , c(1,2,3,4,7,8,9,10,11,5,6)]  
    for(j in 1:length(pil2$POS)){
      read <- pil2$READ[j]
      # Handle forward and reverse reads in the same way 
      read <- gsub('.',',',read, fixed = T)
      read <- gsub('a','A', read)
      read <- gsub('c','C', read)
      read <- gsub('g','G', read)
      read <- gsub('t','T', read)
      # Remove the markers for the start and ends of reads
      read <- gsub('^','',read)
      read <- gsub('$','',read)
      # Remove the deletions/insertions 
      
      # Count the number of bases present for each nucleotide.
      pil2$A_count[j] <- lengths(regmatches(read, gregexpr("A", read)))
      pil2$C_count[j] <- lengths(regmatches(read, gregexpr("C", read)))
      pil2$G_count[j] <- lengths(regmatches(read, gregexpr("G", read)))
      pil2$T_count[j] <- lengths(regmatches(read, gregexpr("T", read)))
      if(pil2$REF[j] == 'A'){pil2$A_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
      if(pil2$REF[j] == 'C'){pil2$C_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
      if(pil2$REF[j] == 'G'){pil2$G_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
      if(pil2$REF[j] == 'T'){pil2$T_count[j] <- lengths(regmatches(read, gregexpr(",", read)))}
      
      # Determine which nucleotide to use as the consensus.
      a <- pil2$A_count[j]
      c <- pil2$C_count[j]
      g <- pil2$G_count[j]
      t <- pil2$T_count[j]
      
      # Determine which nucleotide has the maximum
      max <- max(a,c,g,t)
      
      # Check that the max is above the threshold.
      if(pil2$COVERAGE[j] == 0){max <- -1}else{
        if(max/pil2$COVERAGE[j] < minAF){max <- -1}}
      
      # Assign the consensus. 
      if(pil2$A_count[j] == max){pil2$Consensus[j] <- 'A'}
      if(pil2$C_count[j] == max){pil2$Consensus[j] <- 'C'}
      if(pil2$G_count[j] == max){pil2$Consensus[j] <- 'G'}
      if(pil2$T_count[j] == max){pil2$Consensus[j] <- 'T'}
      if(max == -1){pil2$Consensus[j] <- 'N'} # If there is no clear maximum, then place an N for consensus.
      if(pil2$COVERAGE[j] < minCoverage){pil2$Consensus[j] <- 'N'} # If the coverage is too low, then place an N for consensus.
    }
    # Make the consensus sequence into a data frame.
    col <- c('fasta')
    consensus <- data.frame(matrix(NA, nrow = 2, ncol = length(col)))
    colnames(consensus) <- col
    consensus$fasta[1] <- paste0('>', vsp,'_2')
    consensus$fasta[2] <- paste(t(pil2$Consensus), collapse = '')
    
    # Save the consensus sequence.
    file <- paste0('./03_Consensus/',vsp,'_consensus2.fasta')
    write.table(consensus, file=file, quote=FALSE, sep='\t', col.names = F, row.names = F)
    
    # Save the data frame with nucleotide positions.
    write.csv(pil2, paste0('./02_Variant-Calling/',vsp,'_base-calls.csv'))
  }
  return(ref)
}

# Filters the mutations to call only the major variants.  
major_variants <- function(vsp, minCoverage, minAF, minMQA, minQUAL){
  # vsp <- ref$sample_id[i]
  vcf1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'), comment.char = '#',sep='\t',header=FALSE)
  colnames(vcf1) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENOME')
  # Separate the columns with semicolons to be individual columns 
  # Determine the colnames of the vcf file.
  temp <- strsplit(as.character(data.frame(vcf1$INFO[1])), ";")
  temp <- t(sapply(temp[as.logical(lengths(temp))], function(a) c(a, rep("",max(lengths(temp))-length(a)))))
  col <- as.character(gsub("=.*","",temp))
  # Add the column names to the vcf data frame.
  for(j in 1:length(col)){
    vcf1[col[j]] <- NA
  }
  # Populate the new column from the INFO
  for(j in 1:length(vcf1$POS)){
    temp <- strsplit(as.character(data.frame(vcf1$INFO[j])), ";")
    temp <- data.frame(gsub(".*=","",t(sapply(temp[as.logical(lengths(temp))], function(a) c(a, rep("",max(lengths(temp))-length(a)))))))
    colnames(temp) <- col
    vcf1[j,11:(11+length(col))] <- temp[1,]
  }
  # Save the variant calls as tsv file.
  write.table(vcf1, paste0('./02_Variant-Calling/',vsp,'_allVariants.tsv'), sep='\t',row.names = F, col.names = T)
  vcf1 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.tsv'), sep='\t',header=T)
  vcf2 <- vcf1
  
  # Filter the variant calls by minimum coverage, allele frequence, mapping quality, and base quality 
  vcf2$MQA<- as.numeric(vcf2$MQS)/as.numeric(vcf2$DP)/as.numeric(vcf2$AF) # Get the average mapping quality score.
  vcf2 <- subset(vcf2, as.numeric(vcf2$DP) >= minCoverage)
  vcf2 <- subset(vcf2, as.numeric(vcf2$AF) >= minAF)
  vcf2 <- subset(vcf2, as.numeric(vcf2$MQA) >=minMQ)
  # vcf2 <- subset(vcf2, vcf2$QUAL >= minQUAL) # Commented out because the qual score is not correct.
  write.table(vcf2, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.tsv'), sep='\t',row.names = F, col.names = T)
  vcf2 <- vcf2[,1:10]
  
  # Prepare the headers to save the vcf file
  # vcf3 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'),sep=';',header=FALSE)
  vcf3 <- read.csv2(paste0('./02_Variant-Calling/',vsp,'_allVariants.vcf'), header = F,sep = "!", quote="")
  vcf3 <- subset(vcf3,grepl('#',vcf3$V1))
  temp <- data.frame(rbind(rep(NA, ncol(vcf2)),rep(NA, ncol(vcf2))))
  n <- length(vcf3$V1)
  temp <- as.data.frame(lapply(temp, rep, n))[1:n,]
  colnames(temp) <- colnames(vcf2)
  vcf4 <- rbind(temp,vcf2)
  vcf4$CHROM[1:n] <- vcf3$V1
  vcf4[is.na(vcf4)] <-  ''
  
  # Save the filtered vcf file.
  write.table(vcf4, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), sep = '\t', quote = F, row.names = F, col.names = F)
  # vcf5 <- read.csv(paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), header = F, sep = '!')
  vcf5 <- read.csv2(paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), header = F,sep = "!", quote="")
  vcf5 <- data.frame(gsub('\t\t\t\t\t\t\t\t\t','',vcf5$V1))
  write.table(vcf5, paste0('./02_Variant-Calling/',vsp,'_filteredVariants.vcf'), sep = '', quote = F, row.names = F, col.names = F)
}

# Rename the consensus sequence so that it is the VSP.
rename_consensus <- function(vsp){
  # vsp <- ref$sample_id[i]
  fasta <- read.csv2(paste0('./03_Consensus/',vsp,'_consensus.fasta'), header = F,sep = "!", quote="")
  fasta$V1[1] <- paste0('>',vsp)
  write.table(fasta, paste0('./03_Consensus/',vsp,'_consensus.fasta'), sep = '', quote = F, row.names = F, col.names = F)
}

# Determine the sequencing statistics from the output files and record them to the reference file.  
extract_stat <- function(i, ref){
  # Read in the mapping stat file.
  stat <- read.csv('temp.txt', header = F, sep = '!')
  
  # Get the percent reads mapped
  temp <- subset(stat,grepl('Percent mapped',stat$V1))
  temp <- sub('.*:', '', temp)
  temp <- gsub(' ', '', temp, fixed = TRUE)
  temp <- gsub('\t', '', temp, fixed = TRUE)
  ref$percent_aligned[i] <- as.numeric(temp)
  ref$percent_aligned[is.na(ref$percent_aligned)] <-  0
  
  # Get the average coverage
  temp <- subset(stat,grepl('Average coverage:',stat$V1))
  temp <- sub('.*:', '', temp)
  temp <- gsub(' ', '', temp, fixed = TRUE)
  temp <- gsub('\t', '', temp, fixed = TRUE)
  ref$average_coverage[i] <- as.numeric(temp)
  
  # Get the read 50x coverage for the 550bp target region.
  # Construct the file path
  file_path <- paste0('./02_Variant-Calling/', ref$sample_id[i], '_basecov.txt')
  # Check if the file exists
  if (file.exists(file_path)) {
    # Read the file if it exists
    temp <- read.csv(file_path, sep = '\t')
    temp <- subset(temp, temp$Coverage >= 50)
    nt_length <- 1
    if(ref$amplicon[i]=='RBD'){nt_length <- 550}
    if(ref$amplicon[i]=='3CL-Pro'){nt_length <- 554}
    ref$percent50x_coverage[i] <- round(length(temp$Coverage)/nt_length*100,2)
    print("File has been successfully read.")
  } else {
    # Print a message if the file does not exist
    print("File does not exist.")
  }

  
  return(ref)
}

# Generate the PDF report for each of the samples.
make_report <- function(i, ref){
  
  # Pick the correct start position depending on the amplicon selected
  start_pos <- 0
  if(ref$amplicon[i]=='RBD'){start_pos <- 22772}
  if(ref$amplicon[i]=='3CL-Pro'){start_pos <- 10136}
  
  ## Prepare the third page of the report. This is the first page that is made because the mutation file is required to predict the lineages for the first page.
  # Only prepare the third page if a consensus sequence was made.
  vsp <- ref$sample_id[i]
  
  ## Make the first page of the report.
  sample <- data.frame(t(ref[i,]))
  colnames(sample) <- c(' ')
  
  table <- sample %>%
    kbl(caption = ' ', 
        font_size = 30)  %>%
    kable_classic(full_width = F, html_font = "Cambria") %>%
    kable_styling(latex_options="scale_down")
  table <- add_header_above(table, c('Phage Sample Report' = 2), font_size = 30)
  
  # XXX 
  # Might want to save table as a pdf to be part of the report so that the version number and commands are known.
  
  ## Make the second page of the report.
  file_path <- paste0('./02_Variant-Calling/', ref$sample_id[i], '_basecov.txt')
  
  # Check if the file exists
  if (file.exists(file_path)) {
    
    
    cov1 <- read.csv(paste0('./02_Variant-Calling/',ref$sample_id[i],'_basecov.txt'),sep='\t')
    cov1$Pos <- cov1$Pos + 1
    
    cov2 <- head(cov1,floor(nrow(cov1)/2))
    cov2 <- head(cov1,500000)
    
    # Save the coverage information as a csv file.
    write.csv(cov2,paste0('./04_Reports/table_coverage_',ref$sample_id[i],'.csv'), quote = F, row.names = F)
    
    # Generate the report pdf pages.
    # pdf('temp_page2report.pdf', paper = 'a4', height = 0, width = 0) 
    pdf(paste0('./04_Reports/plot_coverage_',ref$sample_id[i],'.pdf'), height = 17.83, width = 13.78)
    par(mfrow=c(4,1),xpd = T, mar = c(8,8,8,8))
    
    # Prevent scientific notation
    options(scipen = 999)
    
    # Create the barplot
    barplot(cov2$Coverage, border = 'grey40', col = 'grey40', xaxs = 'i', yaxs = 'i', 
            cex.sub = 2, cex.axis = 1.8, cex.lab = 2, cex.names = 1.62, 
            xlab = 'Genome position', ylab = 'Coverage', 
            names.arg = seq(1, nrow(cov2)) + start_pos, xaxt = 'n',
            main = ref$sample_id[i])
    
    dev.off()
    
    
    # Read the file if it exists
    print("File has been successfully read.")
  } else {
    # Print a message if the file does not exist
    print("File does not exist.")
  }
  
  
  return(ref)
}

################################################################################
# Initialization
################################################################################

# Make the necessary file architecture.
directories <- c('00_Raw-Data','temp_Sam-Files','01_Bam-Files','02_Variant-Calling','03_Consensus','04_Reports','05_Summaries')
for(i in 1:length(directories)){
  if(!file.exists(directories[i])){dir.create(directories[i])}
}

# Gather all of the files in the directory.
files <- list.files("./00_Raw-Data/")

# Filter files to only the fastq files. 
files1 <- data.frame(files)
files1$filter <- "n"
for(i in 1:length(files1$files)){
  if(grepl("fastq",files1$files[i]) && grepl("R1",files1$files[i])){
    files1$filter[i] <- "y"
  }
}
files2 <- subset(files1, files1$filter == "y")

# Make a blank data frame
col <- c('PSP', 'sample_id','trial_id','files','sample_type','rationale','sample_date','files','filter')
ref <- data.frame(matrix(NA, nrow = nrow(files2), ncol = length(col)))
colnames(ref) <- col
ref$files <- files2$files
ref$filter <- files2$filter
rownames(ref) <- NULL
ref$files <- gsub('R1','RX',ref$files)

# Save the VSP of the sample as the sample name.
ref$sample_id <- substr(ref$files, 1, nchar(ref$files) - 16)

# # Pull the metadata.
metadata <- read.csv(file_met)
# Determine the sample dates.
ref$sample_date <- NA
ref$trial_id <- NA
ref$rationale <- NA
ref$sample_type <- NA
ref$PSP <- NA
for(ii in 1:length(ref$sample_date)){
  temp <- subset(metadata,metadata$sample_id == ref$sample_id[ii])
  if(nrow(temp)>0){
    ref$PSP[ii] <- as.character(temp$PSP[1])
    ref$sample_date[ii] <- as.character(temp$sample_date[1])
    ref$trial_id[ii] <- as.character(temp$trial_id[1])
    ref$rationale[ii] <- as.character(temp$rationale[1])
    ref$sample_type[ii] <- as.character(temp$sample_type[1])
  }
}
# Determine the report generation date.
ref$report_date <- format(Sys.Date(), "%m/%d/%Y")

# Initialize the reference file.
ref$amplicon <- NA
ref$reads_aligned <- 0
ref$percent_aligned <- 0
ref$average_coverage <- 0
ref$percent50x_coverage <- 0
ref$number_mutations <- NA
ref$predicted_variant <- NA
ref$error <- ''
ref$command_sam <- ''
ref$command_bam <- ''
ref$command_filt1 <- ''
ref$command_qual <- ''
ref$command_sort <- ''
ref$command_index <- ''
ref$command_pileup1 <- ''
ref$command_pileup2 <- ''
ref$command_var1 <- ''
ref$command_var2 <- ''
ref$command_consensus <- ''
ref$version_bwa <- ''
ref$version_samtools <- ''
ref$version_bcftools <- ''
ref$version_bbmap <- ''
# Remove the "filter" column since this is no longer needed.
ref <- ref[ , -which(names(ref) %in% c("filter"))]


# Make an index file for the reference.
ref_files <- c('T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17','T18','T19','T20')
for(ii in 1:length(ref_files)){
  system(paste0('cp ./Bin/',ref_files[ii],'.fasta ./'))
  system(paste0('bwa index ',ref_files[ii],'.fasta'))
}


## Recording Version History.
# BWA
system('bwa 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bwa <- gsub('Version: ','', version$V1[2])
# Samtools
system('samtools 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_samtools <- gsub('Version: ','', version$V1[2])
# Bcftools
system('bcftools 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bcftools <- gsub('Version: ','', version$V1[2])
# BBmap
system('./Bin/bbmap/bbmap.sh -version 2> temp.txt')
version <- read.table('temp.txt', sep = '!')
ref$version_bbmap <- gsub('BBMap version ','', version$V1[2])

# Determine which pipeline to send the samples to.
ref$amplicon <- ref$sample_type

################################################################################
# Execute commands for each sample
################################################################################

# Run the commands for each file.


# Iterate through the contigs within the metadata file.

library(foreach)
library(doParallel)

# Register the number of cores
num_cores_par <- detectCores() - 1
registerDoParallel(cores = num_cores_par)

# for(ii in 1:nrow(ref)){
# foreach(ii = 1:nrow(ref)) %dopar% {
for(ii in 1:nrow(ref)){

  # Pick the database.
  c1_t <- ref$amplicon[ii]
  # List all files in the current directory
  c1 <- list.files(pattern = "\\.fasta$")
  # Subset to only include files that start with the string defined in c1_t
  c1 <- c1[grepl(paste0("^", c1_t), c1)]
  
  # Align the reads to the reference.
  ref$command_sam[ii] <- paste0('bwa mem -M ',c1,' 00_Raw-Data/',gsub("RX", "R1", ref$files[ii]),' 00_Raw-Data/', gsub("RX", "R2", ref$files[ii]), ' >  ./temp_Sam-Files/', ref$sample_id[ii],'_genome.sam')
  system(ref$command_sam[ii])
  
  # Convert sam to bam files. (***CONSIDER changing this to use just 1 core "samtools view -@ 1")
  ref$command_bam[ii] <- paste0('samtools view -S -b -h ./temp_Sam-Files/', ref$sample_id[ii],'_genome.sam',' > ./01_Bam-Files/', ref$sample_id[ii],'_genome.bam')
  system(ref$command_bam[ii])
  
  # Filter reads based on size. # XXX Return to this, it should be used to filter the read lengths so that small reads and long reads are removed but this currently just keeps all reads because of incompatibility of command to do the filtering as it is.
  ref$command_filt1[ii] <- paste0('samtools view -F 8 -h ./01_Bam-Files/', ref$sample_id[ii],'_genome.bam', " | awk 'length($10) > ", as.character(minAmpliconLength), " && length($10) < ", as.character(as.integer(maxAmpliconLength)), " || $1 ~ /^@/' | ", 'samtools view -b -h > ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.bam')
  ref$command_filt1[ii] <- paste0('samtools view -h ./01_Bam-Files/', ref$sample_id[ii],'_genome.bam', " | awk 'length($10) > ", as.character(minAmpliconLength), " && length($10) < ", as.character(as.integer(maxAmpliconLength)), " || $1 ~ /^@/' | ", 'samtools view -b -h > ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.bam')
  ref$command_filt1[ii] <- paste0('samtools view -h ./01_Bam-Files/', ref$sample_id[ii],'_genome.bam > ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.bam')
  system(ref$command_filt1[ii])
  
  # Filter reads by mapping quality.
  ref$command_qual[ii] <- paste0('samtools view -h -q ',minMQ,' -b ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.bam > ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.bam')
  system(ref$command_qual[ii])
  
  # Determine the number of reads that aligned.
  lengthAlignedReadsIDs <- length(system(paste0('samtools view ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.bam | cut -f 1 | uniq'), intern = TRUE))
  
  # Report the number of reads that aligned.
  ref$reads_aligned[ii] <- lengthAlignedReadsIDs * 2
  
  # Sort the sample reads.
  ref$command_sort[ii] <- paste0('samtools sort -o ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.sorted.bam ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.bam')
  system(ref$command_sort[ii])
  
  # Index the sample
  ref$command_index[ii] <- paste0('samtools index ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.sorted.bam')
  system(ref$command_index[ii])
  
  # Make the pileup through bbmap.
  file_pileup <- paste0(ref$PSP[ii],'.txt')
  ref$command_pileup1[ii] <- paste0('./Bin/bbmap/pileup.sh in=./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.sorted.bam', ' ref=',c1,' out=./02_Variant-Calling/',ref$sample_id[ii],'_covstats.txt basecov=./02_Variant-Calling/',ref$sample_id[ii],'_basecov.txt countgc=f overwrite=t 2> ',file_pileup)
  system(ref$command_pileup1[ii])
  
  # Extract the stats from bbmap and record them to the stat file.
  ref <- extract_stat(ii, ref)
  system(paste0('rm -r ',file_pileup))
  
  # Make the pileup through samtools.
  ref$command_pileup2[ii] <- paste0('samtools mpileup ./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.sorted.bam -f ./',c1,' > ./02_Variant-Calling/',ref$sample_id[ii],'_pileup.tsv')
  system(ref$command_pileup2[ii])
  
  # Make a consensus file if there are enough reads.
  if(ref$percent50x_coverage[ii] < QCThreshold){ref$error[ii] <- paste0(ref$error[ii],'|no consensus made - too low coverage|')}
  if(!grepl('too low coverage', ref$error[ii])){ref <- make_consensus(ii, minAF, minCoverage, ref)}
  
  # Call variants.
  ref$command_var1[ii] <- paste0('./Bin/bbmap/callvariants.sh in=./01_Bam-Files/',ref$sample_id[ii],'_genome.filt.qual.sorted.bam ref=',c1,' out=./02_Variant-Calling/',ref$sample_id[ii],'_allVariants.vcf shist=./02_Variant-Calling/',ref$sample_id[ii],'_variantQualityHisto.txt rarity=0 overwrite=t clearfilters')
  system(ref$command_var1[ii])
  
  # # If there are enough reads, then attempt to make a consensus.
  # if(!grepl('no consensus made',ref$error[ii])){ # Check to make sure a concensus should be made.
  #   # Filter variant calls
  #   major_variants(ref$sample_id[ii], minCoverage, minAF, minMQA, minQUAL)
  #   
  #   # Compress the variant calls and index them.
  #   ref$command_var2[ii] <- paste0('bcftools view ./02_Variant-Calling/',ref$sample_id[ii],'_filteredVariants.vcf -Oz -o ./02_Variant-Calling/',ref$sample_id[ii],'_filteredVariants.vcf.gz')
  #   system(ref$command_var2[ii])
  #   ref$command_var2[ii] <- paste0('bcftools index ./02_Variant-Calling/',ref$sample_id[ii],'_filteredVariants.vcf.gz')
  #   system(ref$command_var2[ii])
  #   
  #   # # Get the consensus sequence.
  #   ref$command_consensus[ii] <- paste0('cat ',c1,' | bcftools consensus ./02_Variant-Calling/',ref$sample_id[ii],'_filteredVariants.vcf.gz > ./03_Consensus/',ref$sample_id[ii],'_consensus.fasta')
  #   system(ref$command_consensus[ii])
  #   
  #   # Rename the consensus sequence. 
  #   rename_consensus(ref$sample_id[ii])
  # }
  
  # Generate the sample report.
  ref <- make_report(ii, ref)
  
  # Save the metadata file as it is running.
  ref_temp <- ref[ii,]
  # Check if the file exists
  if (!file.exists(file_ref_out)) {
    # If the file does not exist, write the data to the file
    write.csv(ref_temp, file_ref_out, row.names = FALSE)
  } else {
    # If the file exists, append the data to the existing file
    write.table(ref_temp, file_ref_out, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
}

# Record the reference file.
write.csv(ref, './05_Summaries/reference.csv', row.names = F)

################################################################################
# Cleanup
################################################################################
system('rm -r ./nc_045512.2_*')
system('rm -r temp*')
system('rm -r 05_Summaries/temp_*')



# Function that determines all of the coverage tables then makes a single plot that contains all of the information.
get_cov_plot <- function(){
  gcp1 <- list.files('04_Reports/')
  gcp2 <- gcp1[grep('table_coverage',gcp1)] # Pull only the files that are the tables for coverage
  gcp3 <- gsub('table_coverage_','',gsub('.csv','',gcp2)) # Trim the file name to remove the header table_coverage and the file extension.
  
  tt_types <- unique(ref$sample_type)
  
  # Iterate through each of the phage types.
  for(kk in 1:length(tt_types)){
    
    # Determine which files should be used.
    tt <- tt_types[kk]
    t1 <- subset(ref,ref$sample_type == tt)
    t1_csv <- paste0('table_coverage_',t1$sample_id,'.csv')
    t1_short <- gsub('table_coverage_','',gsub('.csv','',t1_csv)) # Trim the file name to remove the header table_coverage and the file extension.
    gcp2t <- t1_csv[t1_csv %in% gcp2]
    gcp3t <- t1_short[t1_short %in% gcp3]
    
    if(length(gcp3t) > 0){
      # Iterate through the samples to determine rgw maximum
      gcp_max <- 0
      for(jj in 1:length(gcp2t)){
        gcp_file <- gcp2t[jj]
        
        cov2 <- read.csv(paste0('./04_Reports/',gcp_file))
        gcp_max <- max(gcp_max,max(cov2$Coverage))
      }
      
      # Make a pdf that contains plots of all the samples.
      file_plot_coverage_summary <- paste0('./04_Reports/plot_coverage_summary_',tt,'.pdf')
      pdf(file_plot_coverage_summary, height = 17.83, width = 13.78)
      par(mfrow=c(4,1),xpd = T, mar = c(8,8,8,8))
      
      for(jj in 1:length(gcp2t)){
        gcp_file <- gcp2t[jj]
        
        cov2 <- read.csv(paste0('./04_Reports/',gcp_file))
        
        # Determine what would be considered the highest XXXth percentile.
        # Here it is the top 5% that are not 0
        cov3 <- subset(cov2,cov2$Coverage != 0)
        threshold <- quantile(cov3$Coverage, 0.95)
        
        # Prevent scientific notation
        options(scipen = 999)
        
        # Create the barplot
        bar_colors <- ifelse(cov2$Coverage > threshold, 'blue', 'black')
        border_colors <- ifelse(cov2$Coverage > threshold, 'blue', 'grey40')
        
        barplot(cov2$Coverage, border = border_colors, col = bar_colors, xaxs = 'i', yaxs = 'i', 
                cex.sub = 2, cex.axis = 1.8, cex.lab = 2, cex.names = 1.62, 
                xlab = 'Genome position', ylab = 'Coverage', 
                names.arg = seq(1, nrow(cov2)) + min(cov2$Pos), xaxt = 'n',
                main = gcp3t[jj], ylim = c(0, gcp_max))
      }
      dev.off()
    }
  }
}
get_cov_plot()

# Function that checks which files have not been written yet.
dir_report <- '04_Reports'
get_incomplete_samples <- function(ref,dir_report){
  # Determine which files were created.
  file_report <- list.files(dir_report)
  
  # Determine which files were 
  file_report1 <- data.frame(file_report)
  file_report2 <- data.frame(file_report1[grepl(pattern = "table_coverage_", x = file_report1$file_report, ignore.case = TRUE), ])
  colnames(file_report2) <- 'file'
  file_report2$sample <- gsub('.csv','',gsub('table_coverage_','',file_report2$file))
  
  # Determine those samples that are missing.
  ref2 <- subset(ref, !ref$sample_id %in% file_report2$sample)
  
  return(ref2)
}

# ref_original <- ref
ref <- get_incomplete_samples(ref_original,dir_report)

