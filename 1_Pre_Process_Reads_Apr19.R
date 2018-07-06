###### This script removes barcodes, decontaminates and merges overlapping reads
###### Requires bbmap tools and afterQC

##########################################################################################################
#Step 1: Set up everything for the QC run
##########################################################################################################

#Cluster
raw.dir<-"/home/c111h652/scratch/Reduced_Probes/raw_reads"
work.dir<-"/home/c111h652/scratch/Reduced_Probes"

#Directory of the raw reads
after.path<-"/home/c111h652/programs/AfterQC"
adapter.path<-"/home/c111h652/programs/bbmap/resources/adapters.fa"
decontam.path<-"/home/c111h652/scratch/Contamination_Genomes"

#HOME 
#Directory locations
#work.dir<-"/Volumes/Armored/Mantellidae_Dataset"
#raw.dir<-"/Volumes/Armored/Raw_Data/Mantellidae_Dataset"

#Paths to important things used
#after.path<-"/Programs/AfterQC" #After QC path if R isn't recognizing 
#decontam.path<-"/Volumes/Armored/Genetic_Resources/Contamination_Genomes" 
#adapter.path<-"/Programs/bbmap/resources/adapters.fa" #adapter path in bbmap

#Various settings of choice
pypy<-"no" #if installed, speeds up AfterQC
#Sample file
# "File" column for you file name
# "Sample" column should correspond to the file name for each samples set of reads e.g. "CRH01"
sample.file<-"File_rename.csv" #two column file with "File" and "Sample"
out.dir<-"Processed_Samples"
threads<-6
mem<-"80g"

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
  #Finds all files for this given sample and turns into a cat string
  setwd(work.dir)
  sample.reads<-reads[grep(pattern = sample.data$File[i], x = reads)]
  
  if (length(sample.reads) == 0){ sample.reads<-reads[grep(pattern = sample.data$Sample[i], x = reads)] }
  
  if (length(sample.reads) == 0 ){ 
    print(paste(sample.data$Sample[i], " does not have any reads present for files ", 
                sample.data$File[i], " from the input spreadsheet. ", sep = ""))
    next
    } #end if statement
  
  #Create sample directory
  dir.create(sample.data$Sample[i])
 
  #Create symlinks to raw reads
  dir.create(paste(sample.data$Sample[i], "/raw_reads_symlink", sep = ""))
  setwd(paste(sample.data$Sample[i], "/raw_reads_symlink", sep = ""))
  system(paste("ln -s ", raw.dir, "/", sample.reads[1], " ", sample.data$Sample[i], "_READ1.fastq.gz", sep = ""))
  system(paste("ln -s ", raw.dir, "/", sample.reads[2], " ", sample.data$Sample[i], "_READ2.fastq.gz", sep = ""))
  
  dir.create(paste(work.dir, "/", sample.data$Sample[i], "/cleaned_reads_snp", sep = ""))
  setwd(paste(work.dir, "/", sample.data$Sample[i], "/cleaned_reads_snp", sep = ""))
  
  #Runs bbduk, which removes the easy adapter contamination.
  system(paste("bbduk.sh in1=", raw.dir, "/", sample.reads[1], " in2=", raw.dir, "/", sample.reads[2], 
               " ref=", adapter.path, " ftm=5 ktrim=r k=23 mink=8 hdist=1 tbo tpe minlength=25",
               " outs=clean_singletons.fastq.gz out1=clean_READ1.fastq.gz out2=clean_READ2.fastq.gz", sep = ""))
  
   #Next runs bbsplit to remove other sources of contamination from other organisms
  system(paste("bbsplit.sh -Xmx", mem," in1=clean_READ1.fastq.gz in2=clean_READ2.fastq.gz",
               " ref=",decontam.path, " minid=0.95",
               " outu1=", sample.data$Sample[i], "_READ1.fastq.gz",
               " outu2=", sample.data$Sample[i], "_READ2.fastq.gz", sep = ""))
  
  #Decontaminate reads singletons
  system(paste("bbsplit.sh -Xmx", mem," in1=", work.dir, "/", sample.data$Sample[i], 
               "/cleaned_reads_snp/clean_singletons.fastq.gz", 
               " outu1=", sample.data$Sample[i], "_singletons.fastq.gz minid=0.95", sep = ""))
  
  #Deletes the temp files
  system("rm clean_READ1.fastq.gz clean_READ2.fastq.gz clean_singletons.fastq.gz")
  system("rm -r ref")
} #end i loop for first round adapter removal


##########################################################################################################
#Step 3: Use AfterQC to remove adapters / error correct
##########################################################################################################

setwd(work.dir)
sample.reads<-list.files(".", full.names = F, recursive = T)
sample.reads<-sample.reads[grep(pattern = "cleaned_reads_snp", x = sample.reads)]
sample.reads<-sample.reads[grep(pattern = "READ1|READ2", x = sample.reads)]
dir.create("AQC_run/")

for (i in 1:length(sample.reads)) {
  #Finds all files for this given sample and turns into a cat string
  setwd(work.dir)
  
  #Gets the number of reads according to threads
  aqc.reads<-sample.reads[1:(2*floor(threads/2))]
  aqc.reads<-aqc.reads[is.na(aqc.reads) !=T]
  
  #Copy files to temp directory
  system(paste("cp ", paste(aqc.reads, collapse = " "), " AQC_run/", sep = ""))

  if (pypy == "yes" | pypy == "Yes"){
      system(paste("pypy ", after.path, "/after.py --read1_flag READ1 --read2_flag READ2 -d AQC_run/",
               " -g good -b bad -q 0 -n 10 -s 35 -u 60 -f 0 -t 0 --trim_pair_same=false", sep = ""))
  } else {
      system(paste("python ", after.path, "/after.py --read1_flag READ1 --read2_flag READ2 -d AQC_run/",
               " -g good -b bad -q 0 -n 10 -s 35 -u 60 -f 0 -t 0 --trim_pair_same=false", sep = ""))
  } #end pypy if
  
  #Break loop when done
  system(paste("rm -r AQC_run/"))
  dir.create("AQC_run/")
  sample.reads<-sample.reads[!sample.reads %in% aqc.reads]
  if (length(sample.reads) == 0) { 
    system(paste("rm -r AQC_run/"))
    system("rm -r bad")
    break }
}#end AQC Loop

##########################################################################################################
#Step 4: Merge paired end reads and decontamination
##########################################################################################################

#Gets final folders
qc.files<-list.files("QC", full.names = T, recursive = T)
setwd(paste(work.dir, "/good", sep = ""))
sample.reads<-list.files(".", full.names = F, recursive = F)

for (i in 1:nrow(sample.data)) {
  
  #Gathers read data and makes new folder for assembly ready reads
  setwd(work.dir)
  good.reads<-sample.reads[grep(pattern = sample.data$Sample[i], x= sample.reads)]
  
  #Skips if the reads are not present or it matches to more than 1 individual
  if (length(good.reads) != 2){
    print(paste(sample.data$Sample[i], " encountered a problem where READ1 and/or READ2 are absent.", sep = ""))
    next
  } #end if statement
  
  #Moves the QC file to a new QC folder
  sample.qc<-qc.files[grep(pattern = sample.data$Sample[i], x= qc.files)]
  dir.create(paste(work.dir, "/", sample.data$Sample[i], "/QC_files", sep = ""))
  system(paste("cp ", paste(sample.qc, collapse = " "), " ", sample.data$Sample[i], "/QC_files", sep = ""))
  
  #Duplicate removal
  setwd(paste(work.dir, "/good", sep = ""))
  save.dir<-paste(work.dir, "/", sample.data$Sample[i], "/assembly_reads", sep = "")
  dir.create(save.dir)

  #paired end read merging
  system(paste("bbmerge-auto.sh in1=", good.reads[1], " in2=", good.reads[2],
               " verystrict=t rem k=60 extend2=60 ecct", 
               " outu=", save.dir, "/", sample.data$Sample[i], "_READ1.fastq.gz",
               " outu2=", save.dir, "/", sample.data$Sample[i], "_READ2.fastq.gz",
               " out=", save.dir, "/merged.fastq.gz", sep = ""))
  
  #combine the singletons file from filtering before with the merged reads
  system(paste("cat ", work.dir, "/", sample.data$Sample[i], "/cleaned_reads_snp/", sample.data$Sample[i], "_singletons.fastq.gz ",
               save.dir, "/merged.fastq.gz  > ", save.dir,"/", sample.data$Sample[i], "_singletons.fastq.gz", sep = ""))
  #remove merged file
  system(paste("rm ", save.dir, "/merged.fastq.gz", sep = ""))
  
}#End merge and decontam loop

#remove the work folder
setwd(work.dir)
system("rm -r good")
system("rm -r QC")

#Should have a folder structure of:
# Preprocess (working directory)
#   1. -> Sample_Folder
#     a. -> raw_reads_symlink
#     b. -> cleaned_reads_snp
#     c. -> assembly_reads
#     d. -> QC_files



