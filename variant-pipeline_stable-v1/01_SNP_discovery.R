##########################################
### NOTE: Run script 1 from bioinformatics-pipeline
##########################################

library(stringr)
library(data.table)
library(ShortRead)
library(Rsamtools)
library(vcfR)
library(VariantAnnotation)
library(TVTB)
library(seqinr)

#Parallelization
library(parallel)
library(foreach)
library(doParallel)

options(stringsAsFactors = FALSE)
#options(warn=0) #for debugging warnings in loops

#Function that makes consensus sequences 
make.consensus<- function (input.alignment, method = c("majority", "threshold", "IUPAC", 
                                                       "profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = c("DNA", "RNA")) {
  
  #input.alignment<-trimmed
  #Converts alignment to matrix of characters to be used
  new.align<-strsplit(as.character(input.alignment), "")
  align.in<-matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)
  
  #Does based on method
  method <- match.arg(method)
  
  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC, 
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold, 
                  res, NA)
    names(res) <- NULL
  }
  
  out.consensus<-DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}

############################################################################################
########### Setup  #########################################################################
##### Setup for running SNP discovery. 
############################################################################################

###### Directory setup ########
threads = 8 # number of cores
mem.val = 80 #in GB
auto.readgroup = TRUE #You want this done automatically from the fastq files. False just ignores this if it crashes
reference.type = "alignment" #"consensus" = consensus seq from alignments or #"sample" uses the sample's contigs

#Setup directories
work.dir<-"/your/work/directory/SNP"
proc.dir<-"/your/work/directory/Processed_Samples"
align.dir<-"/your/work/directory/Alignments"
run.name<-"Output-Name"

############################################################################################
########### Step 1 #########################################################################
##### Make the reference
############################################################################################

#Go to wd and grab files and get sample names
setwd(proc.dir)
sample.datasets<-list.dirs(proc.dir, recursive = F)

##### Alignment ###############
############################

if (reference.type == "alignment"){
  
  ref.path<-paste0(work.dir, "/Reference")
  if( dir.exists(ref.path) == T) { system(paste0("rm -r ", ref.path)) }
  dir.create(ref.path)
  
  setwd(work.dir)
  
  #Go to wd and grab files and get sample names
  sample.datasets<-list.dirs(proc.dir, recursive = F)
  sample.taxa<-gsub(".*/", "", sample.datasets)
  
  #Gathers alignments
  setwd(paste0(align.dir))
  locus.names<-list.files(".")
  
  #Loops through each locus and does operations on them
  final.con<-c()
  for (i in 1:length(locus.names)){
    #Reads in files
    red.align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))

    if (length(red.align) == 0){ next }
    
    #Get and save consensus sequence
    temp.cons<-DNAStringSet(unlist(as.character(red.align)))
    names(temp.cons)<-names(red.align)
    
    #Makes consensus 
    #con.seq<-ConsensusSequence(temp.cons, ambiguity = F, threshold = 0.95, ignoreNonBases = T, noConsensusChar = "N")
    #Gets consensus seq for trimming more
    con.seq<-make.consensus(temp.cons, method = "majority")
    #Removes the edge gaps
    con.seq<-gsub("\\+", "N", as.character(con.seq))
    
    #Removes the edge gaps
    ref.aligned<-as.character(con.seq)
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    
    if (length(not.gaps) <= 50){ next }
    
    red.seq<-subseq(con.seq, start = as.numeric(not.gaps[1]), end = as.numeric(not.gaps[length(not.gaps)]) )
    red.seq<-DNAStringSet(gsub(pattern = "-", replacement = "", x= as.character(red.seq)))
    
    if (width(red.seq) <= 50){ next }
    
    #Gives save name and saves
    names(red.seq)<-gsub(".phy", "", locus.names[i])
    final.con<-append(final.con, red.seq)
  }#end i loop
  
  #Saves final set
  setwd(paste0(work.dir, "/Reference"))
  final.loci<-as.list(as.character(final.con))
  seqinr::write.fasta(sequences = final.loci, names = names(final.loci), 
                      paste("reference.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  #Indexes the reference
  system(paste0("bwa index -a bwtsw reference.fa"))
  system(paste0("samtools faidx reference.fa"))
  system(paste0("picard -Xmx", mem.val, "G CreateSequenceDictionary",
                " REFERENCE=reference.fa OUTPUT=reference.dict"))
}#  consensus reference end


##### Consensus ###############
############################

if (reference.type == "consensus"){

  ref.path<-paste0(work.dir, "/Reference")
  if( dir.exists(ref.path) == T) { system(paste0("rm -r ", ref.path)) }
  dir.create(ref.path)
  
  setwd(work.dir)
  
  #Go to wd and grab files and get sample names
  sample.datasets<-list.dirs(proc.dir, recursive = F)
  sample.taxa<-gsub(".*/", "", sample.datasets)
  
  #Gathers alignments
  setwd(paste0(align.dir))
  locus.names<-list.files(".")
  
  #Loops through each locus and does operations on them
  final.con<-c()
  for (i in 1:length(locus.names)){
    #Reads in files
    align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))
    red.align<-align[names(align) %in% sample.taxa]
    
    if (length(red.align) == 0){ next }
   
    #Get and save consensus sequence
    temp.cons<-DNAStringSet(unlist(as.character(red.align)))
    names(temp.cons)<-names(red.align)
    
    #Makes consensus 
    #con.seq<-ConsensusSequence(temp.cons, ambiguity = F, threshold = 0.95, ignoreNonBases = T, noConsensusChar = "N")
    #Gets consensus seq for trimming more
    con.seq<-make.consensus(temp.cons, method = "majority")
    #Removes the edge gaps
    con.seq<-gsub("\\+", "N", as.character(con.seq))
    
    #Removes the edge gaps
    ref.aligned<-as.character(con.seq)
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    
    if (length(not.gaps) <= 50){ next }
    
    red.seq<-subseq(con.seq, start = as.numeric(not.gaps[1]), end = as.numeric(not.gaps[length(not.gaps)]) )
    red.seq<-DNAStringSet(gsub(pattern = "-", replacement = "", x= as.character(red.seq)))
    
    if (width(red.seq) <= 50){ next }
    
    #Gives save name and saves
    names(red.seq)<-gsub(".phy", "", locus.names[i])
    final.con<-append(final.con, red.seq)
  }#end i loop
  
  #Saves final set
  setwd(paste0(work.dir, "/Reference"))
  final.loci<-as.list(as.character(final.con))
  seqinr::write.fasta(sequences = final.loci, names = names(final.loci), 
              paste("reference.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  #Indexes the reference
  system(paste0("bwa index -a bwtsw reference.fa"))
  system(paste0("samtools faidx reference.fa"))
  system(paste0("picard -Xmx", mem.val, "G CreateSequenceDictionary",
                " REFERENCE=reference.fa OUTPUT=reference.dict"))
}#  consensus reference end


##### Sample ###############
############################

if (reference.type == "sample"){
  
  for (i in 1:length(sample.datasets)){
    
    #Makes reference for each sample
    ref.path<-paste0(sample.datasets[i], "/Reference")
    if( dir.exists(ref.path) == T) { system(paste0("rm -r ", ref.path)) }
    dir.create(ref.path)
    
    #Reads in files
    setwd(ref.path)
    sample.name<-gsub(".*\\/", "", sample.datasets[i])
    #Reads in the contig file to check
    sample.contigs<-paste0(sample.datasets[i], "/assembled-contigs/", sample.name, "_dipcontigs.fa")
    
    #Copies ortholog file to act as a reference
    sample.name<-gsub(".*/", "", sample.contigs)
    system(paste0("cp ", sample.contigs, " ", ref.path, "/reference.fa"))
    system(paste0("bwa index reference.fa"))
    system(paste0("samtools faidx reference.fa"))
    system(paste0("picard -Xmx", mem.val, "G", 
                  " CreateSequenceDictionary REFERENCE=reference.fa OUTPUT=reference.dict"))
    
  }#end i loop
  
}#end if statement

############################################################################################
########### Step 2 #########################################################################
##### Matches Raw reads to reference contigs and creates BAM file. Cleans BAM file. 
############################################################################################

#Sets up multiprocessing
cl = makeCluster(threads)
registerDoParallel(cl)
mem.cl<-floor(mem.val/threads)

#Loops through each locus and does operations on them
foreach(i=1:length(sample.datasets), .packages = c("foreach", "ShortRead")) %dopar% {
#Loops through each locus and does operations on them
#for (i in 1:length(sample.datasets)){
  
  #Reads in files
  setwd(sample.datasets[i])
  sample.name<-gsub(".*\\/", "", sample.datasets[i])

  #Gets reference path
  if (reference.type == "sample"){ ref.path<-paste0(sample.datasets[i], "/Reference") }
  
  #Gets reads together for this sample
  sample.reads<-list.files("cleaned-reads-snp")
  read1<-paste0(sample.datasets[i], "/cleaned-reads-snp/", sample.reads[grep("READ1", sample.reads)])
  read2<-paste0(sample.datasets[i], "/cleaned-reads-snp/", sample.reads[grep("READ2", sample.reads)])

  #Sets up save place
  dir.create(paste0("snp-analysis-", reference.type))
  setwd(paste0(sample.datasets[i], "/snp-analysis-", reference.type))
  
  #Create unmapped reads set
  ############################
  #convert fastqs to a sam file
  system(paste0("picard -Xmx", mem.cl, "G FastqToSam",
                " FASTQ=", read1, " FASTQ2=", read2, " OUTPUT=fastqsam.bam ",
                " SAMPLE_NAME=", sample.name))
  
  #Revert the sam to a bam file. Cleans and compresses
  system(paste0("picard -Xmx", mem.cl, "G RevertSam",
                " I=fastqsam.bam O=revertsam.bam",
                " SANITIZE=true MAX_DISCARD_FRACTION=0.005",
                " ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS",
                " ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname",
                " RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true",
                " REMOVE_ALIGNMENT_INFORMATION=true"))
  
  #Tries to automatically find the read groups from the fasta headers
  if (auto.readgroup == T){
    #Illumina machines generally
    #@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number | barcode1'+barcode2'>
    all.data = id(readFastq(read2))[1]
    tmp.data = as.character(all.data)
    header.data = unlist(strsplit(tmp.data, ":"))
    header.data = header.data[1:4]
    names(header.data) = c("instrument", "run", "flowcell", "lane")
    
    RGID = paste0(header.data["flowcell"], ".", header.data["lane"])
    RGPU = paste0(header.data["flowcell"], ".", header.data["lane"], ".", sample.name)
    
    #Read groups are assigned 
    system(paste0("picard -Xmx", mem.cl, "G AddOrReplaceReadGroups", 
                  " I=revertsam.bam O=unmapped.bam",
                  " RGSM=", sample.name, " RGPU=", RGPU, " RGID=",RGID,
                  " RGLB=LIB-", sample.name, " RGPL=ILLUMINA"))
  } else {
    #Assign read groups all the same
    system(paste0("picard -Xmx", mem.cl, "G AddOrReplaceReadGroups", 
                  " I=revertsam.bam O=unmapped.bam",
                  " RGSM=", sample.name, " RGPU=FLOWCELL1.LANE1.", sample.name," RGID=FLOWCELL1.LANE1",
                  " RGLB=LIB-", sample.name, " RGPL=ILLUMINA"))
  }#end if 
  
  #Checks read group
  #system(paste0("samtools view -H unmapped.bam | grep '@RG'"))
  
  #Intermediate files are deleted
  system("rm fastqsam.bam revertsam.bam")
  
}#end i loop: create 

stopCluster(cl)


############################################################################################
########### Step 3 #########################################################################
##### Aligns with BWA and converts file
############################################################################################

#Go to wd and grab files and get sample names
setwd(proc.dir)
sample.datasets<-list.dirs(proc.dir, recursive = F)

for (i in 1:length(sample.datasets)){
  
  #Get to sample directory
  setwd(paste0(sample.datasets[i], "/snp-analysis-", reference.type))
  sample.name<-gsub(".*\\/", "", sample.datasets[i])
  
  #Gets reference path
  if (reference.type == "sample"){ ref.path<-paste0(sample.datasets[i], "/Reference") }
  
  #Checks to see if file is there
  if(file.exists("unmapped.bam") == F){ error(paste0("BAM file not created for ", sample.datasets[i])) }
  
  #Run BWA Mem
  system("set -o pipefail")
  
  #Piped verison
  system(paste0("picard -Xmx", mem.val, "G SamToFastq", 
                " I=unmapped.bam FASTQ=/dev/stdout TMP_DIR=", sample.datasets[i], "/snp-analysis-",reference.type,"/tmp",
                " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true |",
                " bwa mem -M -p -t ", threads, " ", ref.path, "/reference.fa /dev/stdin |",
                " picard -Xmx", mem.val, "G MergeBamAlignment", 
                " ALIGNED_BAM=/dev/stdin UNMAPPED_BAM=unmapped.bam",
                " OUTPUT=cleaned_final.bam R=", ref.path, "/reference.fa CREATE_INDEX=true ADD_MATE_CIGAR=true",
                " CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true", 
                " MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS",
                " TMP_DIR=",sample.datasets[i], "/snp-analysis-", reference.type, "/tmp"))
  
  system(paste0("rm unmapped.bam"))
  system(paste0("rm -r tmp"))
  
  print(paste0(sample.datasets[i], " done!!"))
  system(paste0("samtools view -H cleaned_final.bam | grep '@RG'"))
  
}#end i loop
  

############################################################################################
########### Step 4 #########################################################################
##### Uses BAM file and gathers SNP variant raw data for the VCF file per sample
############################################################################################

#Sets up multiprocessing
cl = makeCluster(threads)
registerDoParallel(cl)
mem.cl<-floor(mem.val/threads)

#Go to wd and grab files and get sample names
setwd(proc.dir)
sample.datasets<-list.dirs(proc.dir, recursive = F)

#Loops through each locus and does operations on them
foreach(i=1:length(sample.datasets), .packages = c("foreach")) %dopar% {
#for (i in 1:length(sample.datasets)){

  #Reads in files
  setwd(sample.datasets[i])
  sample.name<-gsub(".*\\/", "", sample.datasets[i])
  #Reads in the contig file to check
  sample.reads<-list.files("assembled-contigs")
  #sample.contigs<-paste0(sample.datasets[i], "/assembled-contigs/", sample.name, "_dipcontigs.fa")
  sample.contigs<-paste0(sample.datasets[i], "/assembled-contigs/", sample.reads[grep("dipcontigs", sample.reads)])
  
  #Gets reference path
  if (reference.type == "sample"){ ref.path<-paste0(sample.datasets[i], "/Reference") }
  
  #Goes to SNP directory
  setwd(paste0(sample.datasets[i], "/snp-analysis-", reference.type))
  
  #Sort by coordinate for input into MarkDuplicates
  system(paste0("picard -Xmx", mem.cl, "G SortSam",
                " INPUT=cleaned_final.bam OUTPUT=cleaned_final_sort.bam",
                " CREATE_INDEX=true SORT_ORDER=coordinate"))
  
  #Marks duplicate reads
  system(paste0("picard -Xmx", mem.cl, "G MarkDuplicates",
                " INPUT=cleaned_final_sort.bam OUTPUT=cleaned_final_md.bam",
                " CREATE_INDEX=true METRICS_FILE=duplicate_metrics.txt"))
  
  #Sorts and stuff
  system(paste0("picard -Xmx", mem.cl, "G SortSam",
                " INPUT=cleaned_final_md.bam OUTPUT=/dev/stdout SORT_ORDER=coordinate |",
                " picard -Xmx", mem.cl, "G SetNmAndUqTags",
                " INPUT=/dev/stdin OUTPUT=snp_ready.bam CREATE_INDEX=true R=",ref.path, "/reference.fa"))
  
  #Delete old files to make more space
  system("rm cleaned_final.* cleaned_final_sort* cleaned_final_md.bam")
  
  #Starts to finally look for Haplotypes! *here
  system(paste0("gatk --java-options '-Xmx", mem.cl, "G' HaplotypeCaller", 
                " -R ", ref.path, "/reference.fa -O snp_output.g.vcf.gz -I snp_ready.bam",
                " -ERC GVCF --max-alternate-alleles 3 -bamout snp_call.bam"))
}# end step 3
  
stopCluster(cl)


############################################################################################
########### Step 4 #########################################################################
##### Gather stats
############################################################################################
# 
# #Sets up filtering rules to get counts
# SNPrule = VcfFixedRules(exprs = list(qual20 = expression(QUAL >= 20),
#                                      SNP = expression(as.character(REF) %in% c("A", "T", "G", "C") &
#                                                         as.character(ALT) %in% c("A", "T", "G", "C"))))
# 
# Indrule = VcfFixedRules(exprs = list(INDEL = expression(Biostrings::width(REF) > Biostrings::width(ALT)) )) 
# 
# #Sets up header
# header.data<-c("Sample","TotalVariants", "TotalSNPs", "TotalIndels", "TotalQual", 
#                "TotalPass", "UniqueLoci", "MeanSNP", "sdSNP", "MinSNP", "MaxSNP")
# #Sets up data to collect
# collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.datasets), ncol = length(header.data)))
# setnames(collect.data, header.data)
# collect.data[, Sample:=as.character(Sample)]
# 
# for (i in 1:length(sample.datasets)){
#     
#   #Reads in files
#   setwd(paste0(sample.datasets[i], "/snp-analysis-", reference.type))
#   sample.name<-gsub(".*\\/", "", sample.datasets[i])
# 
#   #Qual score for SNP missing, if needed?
#   tabixVcf = Rsamtools::TabixFile(file = "snp_output.g.vcf.gz")
#   cvcf = VariantAnnotation::readVcf(file = tabixVcf)
#   evcf = VariantAnnotation::expand(x = cvcf, row.names = TRUE)
# 
#   #Gets the counts for SNPS and indels
#   all.stats = summary(evalSeparately(SNPrule, evcf, enclos = .GlobalEnv))
#   ind.stats = summary(evalSeparately(Indrule, evcf, enclos = .GlobalEnv))
#  
#   #Collects the data  
#   set(collect.data, i = as.integer(i), j = match("Sample", header.data), value = gsub(".*/", "", sample.datasets[i]) )
#   set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = as.numeric(all.stats[1]) )
#   set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = as.numeric(all.stats[3]) )
#   set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = as.numeric(ind.stats[2]) )
#   set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = as.numeric(all.stats[2]) )
#   set(collect.data, i = as.integer(i), j = match("TotalPass", header.data), value = as.numeric(all.stats[4]) )
#   
#   #Per locus stats DOES NOT GET THE FILTERED STATS ONLY ORIGINAL. WANT FILTERED
#   loci.names<-gsub("_|_.*", "", names(evcf))
#   loci.counts<-table(loci.names)
#   set(collect.data, i = as.integer(i), j = match("UniqueLoci", header.data), value = length(loci.counts) )
#   set(collect.data, i = as.integer(i), j = match("MeanSNP", header.data), value = mean(loci.counts) )
#   set(collect.data, i = as.integer(i), j = match("sdSNP", header.data), value = sd(loci.counts) )
#   set(collect.data, i = as.integer(i), j = match("MinSNP", header.data), value = min(loci.counts) )
#   set(collect.data, i = as.integer(i), j = match("MaxSNP", header.data), value = max(loci.counts) )
#   
#   print(paste0(sample.datasets[i], " done!!"))
# }
# 
# setwd(proc.dir)
# write.table(collect.data, file = paste0(run.name,"-", reference.type, "-Prelim_Sample_SNP_Summary.txt"), sep = "\t", row.names = F)

#########################
###### END SCRIPT
#########################


  