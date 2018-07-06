###### This script deals with transcriptomes and also getting the probe set rename
library(ape)
library(seqinr)
library(stringr)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
#library(ips)
#library(DECIPHER)
#setup parallel backend to use many processors
library(foreach)
library(doParallel)

options(stringsAsFactors = FALSE)
#options(warn=1) #for debugging warnings in loops

###############################################################################
###############################################################################
######################Step 1 pull out loci and align  #########################
###############################################################################
###############################################################################

setwd("/Volumes/Armored/Ranoidea_Analysis/Master_Contigs/sequencing_run2")
file.names<-list.files(pattern = ".", full.names = F, recursive = F)
file.path<-paste("Master_Contigs/sequencing_run2/", file.names, sep = "")

#Go to wd and grab files and get sample names
setwd("/Volumes/Armored/Ranoidea_Analysis")
#setwd("/home/c111h652/scratch/alignments_ranoidea/Alignment_cerato/aligned")
file.conv<-read.csv("Final_Datasets_2017/contigs_to_raw.csv")

### Data to be collected
parameters<-c("Species", "numberLoci", "meanLength","medianLength","minLength","maxLength",
              "meanCov", "medianCov", "minCov", "maxCov", "sdCov")

species.data<-data.table(matrix(0, nrow = length(file.names), ncol = length(parameters)))
setnames(species.data, parameters)
species.data[, Species:=as.character(Species)]

#Locus data
locus.data<-vector("list", length(file.names))

#Loops through each locus and writes each species to end of file
for (i in 1:length(file.names)){
  
  setwd("/Volumes/Armored/Ranoidea_Analysis")
  
  sample<-file.conv[file.conv$Sample == gsub(pattern = "_final.fa", replacement = "", x = file.names[i]),]
  read1<-paste("/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/cleaned_reads_BB/", sample$Sample, 
               "/split-adapter-quality-trimmed/", sample$Sample, "-READ1.fastq.gz", sep ="")
  read2<-paste("/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/cleaned_reads_BB/", sample$Sample, 
               "/split-adapter-quality-trimmed/", sample$Sample, "-READ2.fastq.gz", sep ="")
  singleton<-paste("/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/cleaned_reads_BB/", sample$Sample, 
                   "/split-adapter-quality-trimmed/", sample$Sample, "-READ-singleton.fastq.gz", sep ="")
  
  #Gets rid of duplicates  
  system(paste("awk '/^>/{f=!d[$1];d[$1]=1}f' ", file.path[i], " > ", file.path[i], "sta", sep = ""))
  
  setwd("/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/seq_depth")
  
  #Maps reads and makes reference
  system(paste("bbwrap.sh -Xmx8g in1=", read1, ",", singleton, " in2=", read2, ",null",   
               " out=", sample$Sample, ".sam append ref=../../",file.path[i], "sta", sep = ""))
  
  #system(paste("bbmap.sh -Xmx8g in1=", read1, " in2=", read2,  
  #             " out=", sample$Sample, ".sam ref=../",file.path[i], "sta", sep = ""))
    
  #Gets coverage statistics for contigs
  system(paste("pileup.sh -Xmx8g in=", sample$Sample, ".sam bincov=", sample$Sample, "_bin.txt binsize=1",
               " covstats=", sample$Sample, "_cov.txt overwrite=true", sep = "")) 
               
  #Reads in coverage file
  cov.table<-fread(file = paste(sample$Sample, "_cov.txt", sep = ""),  sep = "\t", stringsAsFactors = FALSE)   # loads up fasta file
  setnames(cov.table, old = "#ID", new = "cdsName")
  bin.table<-fread(file = paste(sample$Sample, "_bin.txt", sep = ""),  sep = "\t", stringsAsFactors = FALSE)   # loads up fasta file
  
  #Collect data on each locus by aggregating it all
  #Number of loci over 5X, 10X, 20X, and 50x
  save.table<-cbind(Species = gsub(pattern = "_final.fa", replacement = "", x = file.names[i]), cov.table)
  save.table$cdsName<-gsub(pattern = ".* \\|", replacement = "", x = save.table$cdsName)
  locus.data[[i]]<-save.table
  
  #Collect species data
  set(species.data, i = as.integer(i), j = match("Species", parameters), value = gsub(pattern = "_final.fa", replacement = "", file.names[i]))
  set(species.data, i = as.integer(i), j = match("numberLoci", parameters), value = nrow(cov.table))
  set(species.data, i = as.integer(i), j = match("meanLength", parameters), value = mean(cov.table$Length))
  set(species.data, i = as.integer(i), j = match("medianLength", parameters), value = median(cov.table$Length))
  set(species.data, i = as.integer(i), j = match("minLength", parameters), value = min(cov.table$Length))
  set(species.data, i = as.integer(i), j = match("maxLength", parameters), value = max(cov.table$Length))
  
  set(species.data, i = as.integer(i), j = match("meanCov", parameters), value = mean(bin.table$Cov))
  set(species.data, i = as.integer(i), j = match("medianCov", parameters), value = median(bin.table$Cov))
  set(species.data, i = as.integer(i), j = match("minCov", parameters), value = min(bin.table$Cov))
  set(species.data, i = as.integer(i), j = match("maxCov", parameters), value = max(bin.table$Cov))
  set(species.data, i = as.integer(i), j = match("sdCov", parameters), value = sd(bin.table$Cov))
  
  #Delete crap files 
  system(paste("rm /Volumes/Armored/Ranoidea_Analysis/", file.path[i], "sta", sep = ""))
  system("rm -r ref")
  system(paste("rm ", sample$Sample, "*", sep = ""))
  print(species.data[i,])
  print(paste(sample$Sample, " complete!"))
}  

write.data<-species.data[species.data$Species != 0,]
write.table(species.data, file = "/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/seq_depth/Species_data.txt", sep = "\t", row.names = F)

write.locus.data<-rbindlist(locus.data, use.names = T, fill = T, idcol = NULL)
write.table(write.locus.data, file = "/Volumes/Armored/Ranoidea_Analysis/Ranoidea_run2/seq_depth/raw_locus_data.txt", sep = "\t", row.names = F)

  
  
  
  
  
  ##############
  #Step 1: Read in each alignment to get stats for each
  #############
  #get saved names
  #Alignment with introns
  alignment<-scanFa(FaFile(paste("Alignment_cerato/aligned/", file.names[i], sep = "")))   # loads up fasta file
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  new.align<-strsplit(as.character(alignment), "")
  new.align<-lapply(new.align, tolower)
  aligned.set<-as.matrix(as.DNAbin(new.align))
  
  #Alignment with only exons
  alignment<-scanFa(FaFile(paste("Alignment_cerato/no_intron/", file.names[i], sep = "")))   # loads up fasta file
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  new.align<-strsplit(as.character(alignment), "")
  new.align<-lapply(new.align, tolower)
  aligned.exon<-as.matrix(as.DNAbin(new.align))
  
  #Gets taxonomy data
  loci.tax<-taxonomy[taxonomy$Sample_Name %in% rownames(aligned.exon),]
  
  ##############
  #Step 2: Collect species stats! 
  ##############
  set(species.data, i = as.integer(i), j = match("cdsName", species.columns), value = gsub(pattern = ".fa", replacement = "", file.names[i]))

  temp<-ifelse(species.columns %in% loci.tax$Sample_Name, 1, 0)[-1] # this does the trick
  for (k in 1:length(temp)){
  
    set(species.data, i = as.integer(i), j = match(species.columns[k+1], species.columns), temp[k])
  }
  
  ##############
  #Step 3: Collect locus stats! 
  ##############
  
  #Saves data for each locus
  set(final.data, i = as.integer(i), j = match("cdsName", parameters), value = gsub(pattern = ".fa", replacement = "", file.names[i]))
  set(final.data, i = as.integer(i), j = match("isParalog", parameters), value = 0)
  set(final.data, i = as.integer(i), j = match("noSpecies", parameters), value = nrow(aligned.set))
  set(final.data, i = as.integer(i), j = match("noClades", parameters), value = length(unique(loci.tax$Super_Family)))
  set(final.data, i = as.integer(i), j = match("noSalam", parameters), value = length(loci.tax$Super_Family[loci.tax$Super_Family == "Salamandra"]))
  set(final.data, i = as.integer(i), j = match("noBasal", parameters), value = length(loci.tax$Super_Family[loci.tax$Super_Family == "Basal"]))
  set(final.data, i = as.integer(i), j = match("noHyloidea", parameters), value = length(loci.tax$Super_Family[loci.tax$Super_Family == "Hyloidea"]))
  set(final.data, i = as.integer(i), j = match("noRanoidea", parameters), value =length(loci.tax$Super_Family[loci.tax$Super_Family == "Ranoidea"]))
  set(final.data, i = as.integer(i), j = match("noFamilies", parameters), value = length(unique(loci.tax$Family)))
  set(final.data, i = as.integer(i), j = match("noTrans", parameters), value = length(loci.tax$Type[loci.tax$Type == "Transcriptomes"]))
  set(final.data, i = as.integer(i), j = match("noSampled", parameters), value = length(loci.tax$Type[loci.tax$Type == "Sequencing"]))
  
  set(final.data, i = as.integer(i), j = match("alignLength", parameters), value = ncol(aligned.set))
  set(final.data, i = as.integer(i), j = match("mDistance", parameters), value = mean(dist.dna(aligned.set, model = "raw", pairwise.deletion = T), na.rm=T))
  set(final.data, i = as.integer(i), j = match("PISabs", parameters), value = new.pis(aligned.set, abs = T))
  set(final.data, i = as.integer(i), j = match("PISper", parameters), value = new.pis(aligned.set, abs = F))
  
  #Saves data for each locus no intron
  set(final.data, i = as.integer(i), j = match("NI_alignLength", parameters), value = ncol(aligned.exon))
  set(final.data, i = as.integer(i), j = match("NI_mDistance", parameters), value = mean(dist.dna(aligned.exon, model = "raw", pairwise.deletion = T), na.rm=T))
  set(final.data, i = as.integer(i), j = match("NI_PISabs", parameters), value = new.pis(aligned.exon, abs = T))
  set(final.data, i = as.integer(i), j = match("NI_PISper", parameters), value = new.pis(aligned.exon, abs = F))
  #Counter
  if (count == 1000){ 
    print(paste("iteration", i, "complete", Sys.time(), sep = " "))
    count<-0
  } else {
    count<-count+1
  }#end else
}# end loop
#stop cluster
#stopCluster(cl)  

write.table(final.data, file = "Alignment_cerato/aligned_loci_stats_all.txt", row.names = F, sep = "\t")
write.table(species.data, file = "Alignment_cerato/aligned_species_stats_all.txt", row.names = F, sep = "\t")
