###### This script assembles cleaned reads using dipspades

#Need to only add an install of SPADES. Requires older version of SPADES (attached in folder here)

#Should have a folder structure of from Step 1: Preprocess reads
# Preprocess (working directory)
#   1. -> Sample_Folder
#     a. -> raw_reads_symlink
#     b. -> cleaned_reads_snp
#     c. -> assembly_reads
#     d. -> QC_files


##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#Set up parameters
threads<-"8"
mem<-"80"

#Directory where the "Processed_Samples" Folder is located
work.dir<-"/Your/Project/Directory"

########################
#Choose which directory to use for assembly:
#aseembly_reads = recommended assembly reads directory (default)
#cleaned_reads_snp = if the assembly_reads are too aggressively filtered. 
read.dir<-"assembly-reads"
out.dir<-"Assembled_Contigs"
in.dir<-"Processed_Samples"
clean<-"True" #Deletes intermediate files


############################################################################################
########### Step 1 #########################################################################
##### Runs assembly
############################################################################################

#Creates the important working directories
dir.create(paste(work.dir, "/", out.dir, sep = ""))
setwd(paste(work.dir, "/", in.dir, sep = ""))
sample.dirs<-list.dirs(path= ".", full.names = F, recursive = F)

#Cycles through each read pair and does things to them
for (i in 1:length(sample.dirs)) {
  
  #Creates assembly reads folder if not present
  save.assem<-paste(work.dir, "/", in.dir,"/", sample.dirs[i],"/assembled-contigs", sep = "")
  dir.create(save.assem)
  
  #Finds all files for this given sample and turns into a cat string
  setwd(paste(work.dir, "/", out.dir, sep = ""))
  
  #Get location of raw reads per sample
  sample.loc<-paste(work.dir,"/",in.dir,"/",sample.dirs[i],"/",read.dir,"/", sep = "")
  
  #Run SPADES on sample
  system(paste("dipspades.py --pe1-1 ", sample.loc, sample.dirs[i], "_READ1.fastq.gz", 
               " --pe1-2 ", sample.loc, sample.dirs[i], "_READ2.fastq.gz",
               " --merged ", sample.loc, sample.dirs[i], "_singletons.fastq.gz",
               " -o ", sample.dirs[i], " -k 21,33,55,77,99,127 --careful -t ", threads, " -m ", mem, 
               " --hap-assembly --expect-gaps", sep = "")) 
  
  #Dipspades run 
  setwd(paste(work.dir, "/", out.dir, "/", sample.dirs[i], "/dipspades/", sep = "")) 
  #Copies and renames dipcontigs
  system(paste("cp consensus_contigs.fasta ", work.dir, "/", out.dir, sep = ""))
  system(paste("mv ", work.dir, "/", out.dir, "/consensus_contigs.fasta ", 
               work.dir, "/", out.dir, "/", sample.dirs[i], ".fa", sep = ""))
  #Copy to species folder
  system(paste("cp consensus_contigs.fasta ", save.assem, sep = ""))
  system(paste("mv ", save.assem, "/consensus_contigs.fasta ", 
               save.assem, "/", sample.dirs[i], "_dipcontigs.fa", sep = ""))
  
  #Dipspades run 
  setwd(paste(work.dir, "/", out.dir, "/", sample.dirs[i], "/spades/", sep = "")) 

  #Copy to species folder
  system(paste("cp contigs.fasta ", save.assem, sep = ""))
  system(paste("mv ", save.assem, "/contigs.fasta ", 
               save.assem, "/", sample.dirs[i], "_contigs.fa", sep = ""))
  
  #copies to contigs folder
  setwd(paste(work.dir, "/", out.dir, "/", sample.dirs[i], "/spades", sep = ""))
  system(paste("rm -r corrected K21 K33 K55 K77 K99 K127 misc tmp", sep = "")) 
  system(paste("cp -r ", work.dir, "/", out.dir, "/", sample.dirs[i], "/spades ", 
               work.dir, "/", in.dir,"/", sample.dirs[i], sep = ""))

  #Moves remainder to nice place
  setwd(paste(work.dir, "/", out.dir, sep = ""))
  system(paste0("rm -r ", sample.dirs[i]))
  
}#end i loop






