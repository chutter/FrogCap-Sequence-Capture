###### This script removes barcodes, decontaminates and merges overlapping reads
#####################################################################################

#Best, quick way to get everything installed and working:

#These commands need to be entered in the terminal command line, on a computer cluster or whereever. 
#Remove # signs, as they comment out the steps here so R doesn't run them.

#--------------------

#Step 1: Download these two programs and install. Make sure they are in your $PATH

#Fastp: https://github.com/OpenGene/fastp
#bbtools: https://jgi.doe.gov/data-and-tools/bbtools/

#Step 2: Make configuration file. This is the "File_rename.csv" file, which has two columns: File and Sample.
#Column File = the unique string that is part of the file name for the two read pairs. 
#Column Sample = what you want the new name to be called. 

#Step 3: put Contamination_Genomes folder somewhere and note the directory.

#Step 4: Install anaconda, or install a new working environment

# Create conda environment, swap out "c111h652" for your own username. Could vary based on cluster.
### This step creates a new environment wherever  you want

#conda create --prefix $WORK/conda/frogcap
#source activate /panfs/pfs.local/work/bi/c111h652/conda/frogcap

#Gather other dependencies from R 
#=setup channels (has to be in this order)
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --add channels defaults

#conda install -c bioconda fastp
#conda install -c r r-base

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#Cluster
work.dir<-"/home/Project_Name"
raw.dir<-"/home/Project_Name/raw_reads"
decontam.path<-"/home/Project_Name/FrogCap_Files/Contamination_Genomes" #folder of contaminants I gave

#Various settings of choice
#Sample file
# "File" column for you file name
# "Sample" column should correspond to the file name for each samples set of reads e.g. "CRH01"
sample.file<-"File_rename.csv" #two column file with "File" and "Sample"
out.dir<-"Processed_Samples" #I would keep this name. 
threads<-10 #number of threads, 8-10 is a recommended amount
mem<-"80g" #need the "g' character here

##########################################################################################################
#Step 1: Set up directories (do not edit further unless you know what is happening)
##########################################################################################################

#Sets directory and reads in sample file
setwd(work.dir)
sample.data<-read.csv(sample.file, stringsAsFactors = F)

#Checks for missing files if completed
dirs<-list.dirs(path = ".", full.names = F, recursive = F)
done.dirs<-dirs[grep("Assembled_Contigs", dirs)]

if (length(done.dirs) >= 1) {
  con.names<-list.files(pattern = "", full.names = F, recursive = T)
  con.names<-con.names[grep("Assembled_Contigs", con.names)]
  con.names<-gsub("\\.fa", "", con.names)
  con.names<-gsub(".*\\/", "", con.names)
  sample.data<-sample.data[is.na(pmatch(sample.data$Sample, con.names)),]
} #end if

### Get raw reads together
setwd(raw.dir)
reads<-list.files(".", full.names = F, recursive = T)
reads<-reads[grep("md5", reads, invert = T)]

#Checks for undetermined reads and excludes these
if (length(grep(pattern = "Undetermined", x = reads)) != 0){ 
  reads<-reads[reads != reads[grep(pattern = "Undetermined", x = reads)]] }

#Creates output directory
dir.create(paste(work.dir, "/", out.dir, sep = ""))
work.dir<-paste(work.dir, "/", out.dir, sep = "")

##########################################################################################################
#Step 2: Remove adapters and quality checking; decontaminate reads
##########################################################################################################

for (i in 1:nrow(sample.data)) {
  #################################################
  ### Part A: prepare for loading and checks
  #################################################
  #Finds all files for this given sample and turns into a cat string
  setwd(work.dir)
  sample.reads<-reads[grep(pattern = sample.data$File[i], x = reads)]
  
  if (length(sample.reads) == 0){ sample.reads<-reads[grep(pattern = sample.data$Sample[i], x = reads)] }
  
  if (length(sample.reads) == 0 ){ 
    print(paste(sample.data$Sample[i], " does not have any reads present for files ", 
                sample.data$File[i], " from the input spreadsheet. ", sep = ""))
    next
    } #end if statement
 
  #################################################
  ### Part B: Create directories and move files
  #################################################
  #Create sample directory
  dir.create(sample.data$Sample[i])
  dir.create(paste0(work.dir, "/", sample.data$Sample[i], "/reports"))
  #Create symlinks to raw reads
  dir.create(paste0(sample.data$Sample[i], "/raw-reads-symlink"))
  setwd(paste0(sample.data$Sample[i], "/raw-reads-symlink"))
  system(paste0("ln -s ", raw.dir, "/", sample.reads[1], " ", sample.data$Sample[i], "_READ1.fastq.gz"))
  system(paste0("ln -s ", raw.dir, "/", sample.reads[2], " ", sample.data$Sample[i], "_READ2.fastq.gz"))
  raw.reads<-list.files(".")
  raw.path<-paste0(work.dir, "/", sample.data$Sample[i], "/raw-reads-symlink/", raw.reads)
  
  #Gathers stats on initial data
  reads.removed<-data.frame(Source = as.character(), nPairs = as.numeric(), remPairs = as.numeric())
  curr.reads<-as.numeric(system(paste0("zcat < ", raw.path[1], " | echo $((`wc -l`/4))"), intern = T))
  temp.remove<-data.frame(Source = "initial raw reads", nPairs = curr.reads, remPairs = as.numeric(0))
  reads.removed<-rbind(reads.removed, temp.remove)
  
  #################################################
  ### Part C: Creates cleaned reads for SNPS
  #################################################
  #Creates directory
  dir.create(paste0(work.dir, "/", sample.data$Sample[i], "/cleaned-reads-snp"))
  setwd(paste0(work.dir, "/", sample.data$Sample[i], "/cleaned-reads-snp"))
  
  #Runs fastp: only does adapter trimming, no quality stuff
  system(paste0("fastp --in1 ",raw.path[1], " --in2 ", raw.path[2], 
                " --out1 ", raw.reads[1], " --out2 ", raw.reads[2],
                " --length_required 30 --low_complexity_filter --complexity_threshold 30",
                " --html adapter_trimming_fastp.html --json adapter_trimming_fastp.json",
                " --report_title ", sample.data$Sample[i]," --thread ", threads))
  
  system(paste0("cp adapter_trimming* ../reports"))
  system(paste0("rm adapter_trimming*"))
  
  #Saves removal data
  curr.reads<-as.numeric(system(paste0("zcat < ", raw.reads[1], " | echo $((`wc -l`/4))"), intern = T))
  temp.remove<-data.frame(Source = "fastp adapter-complexity filter", nPairs = curr.reads, remPairs = reads.removed[1,2]-curr.reads)
  reads.removed<-rbind(reads.removed, temp.remove)
  
  #Next runs bbsplit to remove other sources of contamination from other organisms
  system(paste0("bbsplit.sh -Xmx", mem," in1=", raw.reads[1], " in2=", raw.reads[2],
               " ref=",decontam.path, " minid=0.95",
               " outu1=Cleaned_READ1.fastq.gz outu2=Cleaned_READ2.fastq.gz"), ignore.stderr = F)
  
  system(paste0("rm -r ref"))
  system(paste0("rm ", raw.reads[1], " ", raw.reads[2]))
  system(paste0("mv Cleaned_READ1.fastq.gz ", raw.reads[1]))
  system(paste0("mv Cleaned_READ2.fastq.gz ", raw.reads[2]))
  
  #Saves removal data
  curr.reads<-as.numeric(system(paste0("zcat < ", raw.reads[1], " | echo $((`wc -l`/4))"), intern = T))
  temp.remove<-data.frame(Source = "decontamination", nPairs = curr.reads, remPairs = reads.removed[2,2]-curr.reads)
  reads.removed<-rbind(reads.removed, temp.remove)
  
  #Totals up all removed
  temp.remove<-data.frame(Source = "total reads removed", nPairs = reads.removed[nrow(reads.removed),2], 
                          remPairs = sum(reads.removed[,3]))
  cleaned.removed<-rbind(reads.removed, temp.remove) 
  
  write.csv(cleaned.removed, file = "../reports/cleaned-reads-snp_stats.csv", row.names = F)
  
  #################################################
  ### Part D: Paired-end read merging and duplicate removal 
  #################################################
  
  save.dir<-paste0(work.dir, "/", sample.data$Sample[i], "/assembly-reads")
  dir.create(save.dir)
  
  # DEDUPE almost exact duplicate removal
  system(paste0("dedupe.sh -Xmx", mem," in1=", raw.reads[1], " in2=", raw.reads[2], " ordered=t overwrite=true",
                " out=", save.dir, "/deduped_reads.fastq.gz minidentity=100"), ignore.stderr = T)
  setwd(save.dir)
  system(paste0("reformat.sh in=deduped_reads.fastq.gz out1=deduped_READ1.fastq.gz out2=deduped_READ2.fastq.gz"), ignore.stdout = T)
  
  #Saves removal data (CURRR HERE)
  curr.reads<-as.numeric(system("zcat < deduped_READ1.fastq.gz | echo $((`wc -l`/4))", intern = T))
  temp.remove<-data.frame(Source = "duplicate removal", nPairs = curr.reads, 
                          remPairs = reads.removed[nrow(reads.removed),2]-curr.reads)
  reads.removed<-rbind(reads.removed, temp.remove)
  
  #paired end read merging
  system(paste0("bbmerge-auto.sh -Xmx", mem, " in1=deduped_READ1.fastq.gz in2=deduped_READ2.fastq.gz",
               " verystrict=t rem k=60 extend2=60 ecct", 
               " outu=", sample.data$Sample[i], "_READ1.fastq.gz",
               " outu2=", sample.data$Sample[i], "_READ2.fastq.gz",
               " out=merged_pe.fastq.gz"))

  #Dedupe the merged reads  
  system(paste0("dedupe.sh -Xmx", mem," in=merged_pe.fastq.gz ordered=t overwrite=true",
                " out=", sample.data$Sample[i], "_singletons.fastq.gz minidentity=98"), ignore.stderr = T)
  
  #Totals up all removed
  temp.remove<-data.frame(Source = "total reads removed", nPairs = reads.removed[nrow(reads.removed),2], 
                          remPairs = sum(reads.removed[,3]))
  reads.removed<-rbind(reads.removed, temp.remove) 
  
  write.csv(reads.removed, file = "../reports/assembly-reads_stats.csv", row.names = F)
  
  #Count the number of reads merged (FIX AND CHECK HERE)
  curr.reads<-as.numeric(system(paste0("zcat < ", raw.reads[1], " | echo $((`wc -l`/4))"), intern = T))
  merge.remove<-data.frame(Source = "merged paired reads", nPairs = curr.reads, 
                          mergePairs = reads.removed[nrow(reads.removed),2]-curr.reads)

  write.csv(merge.remove, file = "../reports/read_merging_stats.csv", row.names = F)
  
  #remove merged file
  system(paste0("rm merged_pe.fastq.gz deduped_*"))
  print(paste0(sample.data$Sample[i], " Complete!"))
} #end i loop 



#Should have a folder structure of:
# Preprocess (working directory)
#   1. -> Sample_Folder
#     a. -> raw_reads_symlink
#     b. -> cleaned_reads_snp
#     c. -> assembly_reads
#     d. -> QC_files



