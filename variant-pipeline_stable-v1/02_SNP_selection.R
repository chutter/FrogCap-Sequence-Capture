library(ape)
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

############################################################################################
########### Step 1 #########################################################################
##### Setup for running SNP discovery. 
############################################################################################


###### parameter setup ########
threads<-8 # number of cores
mem.val<-80 #in GB


#Setup directories
db.name = "SNP_Database"
ref.type = "alignment"

#Setup directories
out.dir<-"/your/work/directory/SNP-Analysis"
ref.dir = "/your/work/directory/Reference"

#From previous script
work.dir<-"/your/work/directory"
proc.dir<-"/your/work/directory/Processed_Samples"
align.dir<-"/your/work/directory/Alignments"
run.name<-"Output-Name"


############################################################################################
########### Step 2 #########################################################################
##### Matches Raw reads to reference contigs and creates BAM file. Cleans BAM file. 
############################################################################################

#Get same vcf files together
setwd(out.dir)
dir.create(db.name)
setwd(paste0(out.dir, "/", db.name))

#Get loci names I guess
probe.loci<-scanFa(FaFile(paste0(ref.dir, "/reference.fa")))
loci.names<-names(probe.loci)

#Get multifile databases together
sample.datasets<-list.dirs(proc.dir, recursive = F)
sample.datasets<-paste0("-V ", sample.datasets, "/snp-analysis-", ref.type, "/snp_output.g.vcf.gz")
vcf.files<-paste0(sample.datasets, collapse = " ")

#Sets up multiprocessing
cl = makeCluster(threads)
registerDoParallel(cl)
mem.val<-floor(mem.val/threads)

#Loops through each locus and does operations on them
foreach(i=1:length(loci.names), .packages = c("foreach")) %dopar% {

#Loops through each locus and does operations on them
#for (i in 1:length(loci.names)){
  
  #Combine them into a single database
  if (dir.exists(loci.names[i]) == T) { system(paste0("rm -r ", loci.names[i])) }
  system(paste0("gatk --java-options '-Xmx", mem.val, "G' GenomicsDBImport ", vcf.files,
                " --genomicsdb-workspace-path ", loci.names[i], " --intervals ", loci.names[i]))
  
  #Combine them into a single database
  system(paste0("gatk --java-options '-Xmx", mem.val, "G' GenotypeGVCFs",
                " -R ", ref.dir, "/reference.fa -V gendb://", loci.names[i],
                " --use-new-qual-calculator true -O ", loci.names[i], ".vcf"))
  
  system(paste0("rm -r ", loci.names[i]))
  
}

stopCluster(cl)


############################################################################################
########### Step 2 #########################################################################
##### Matches Raw reads to reference contigs and creates BAM file. Cleans BAM file. 
############################################################################################

#Sets up filtering rules to get counts
SNPrule = VcfFixedRules(exprs = list(qual20 = expression(QUAL >= 20),
                                     SNP = expression(as.character(REF) %in% c("A", "T", "G", "C") &
                                                        as.character(ALT) %in% c("A", "T", "G", "C"))))

Indrule = VcfFixedRules(exprs = list(INDEL = expression(Biostrings::width(REF) > Biostrings::width(ALT)) )) 


#Sets up header
header.data<-c("Locus", "TotalVariants", "TotalSNPs", "TotalIndels", "TotalQual", "TotalPass")
#Sets up data to collect
collect.data<-data.table(matrix(as.numeric(0), nrow = length(loci.names), ncol = length(header.data)))
setnames(collect.data, header.data)
collect.data[, Locus:=as.character(Locus)]

for (i in 1:length(loci.names)){
  
  #Qual score for SNP missing, if needed?
  cvcf = VariantAnnotation::readVcf(file = paste0(loci.names[i], ".vcf"))
  
  #Catches empty SNPs vcf
  if (dim(cvcf)[1] == 0){
    #Collects the data  
    set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = loci.names[i] )
    set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = 0 )
    set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = 0 )
    set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = 0 )
    set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = 0 )
    set(collect.data, i = as.integer(i), j = match("TotalPass", header.data), value = 0 )
    system(paste0("rm -r ", loci.names[i]))
    
    next
  }
  
  evcf = VariantAnnotation::expand(x = cvcf, row.names = TRUE)
  
  #Gets the counts for SNPS and indels
  all.stats = summary(evalSeparately(SNPrule, evcf, enclos = .GlobalEnv))
  ind.stats = summary(evalSeparately(Indrule, evcf, enclos = .GlobalEnv))
  
  #Collects the data  
  set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = loci.names[i] )
  set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = as.numeric(all.stats[1]) )
  set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = as.numeric(all.stats[3]) )
  set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = as.numeric(ind.stats[2]) )
  set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = as.numeric(all.stats[2]) )
  set(collect.data, i = as.integer(i), j = match("TotalPass", header.data), value = as.numeric(all.stats[4]) )
  
  print(paste0(loci.names[i], " done!!"))

}

write.table(collect.data, file = "Locus_SNP_Summary.txt", sep = "\t", row.names = F)



#This to try instead to get the SNPs or whatever
# gatk SelectVariants \
# -R data/ref/ref.fasta \
# -V gendb://my_database \
# -O combined.g.vcf
# 
# #Can also try this to refine a bit more
# gatk VariantFiltration \
# -V trio.vcf \
# -O trio_VF.vcf \
# --genotype-filter-expression "isHet == 1" \
# --genotype-filter-name "isHetFilter"
# 

#After finding best variants, do base recalibration? 

#########################
###### END SCRIPT
#########################


  