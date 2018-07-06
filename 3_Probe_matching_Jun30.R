#### PREREQUISITES OUTSIDE PROGRAMS ######

#PBLAT: http://icebert.github.io/pblat/
#MAFFT 7.312: https://mafft.cbrc.jp/alignment/software/

#### BIOCONDUCTOR INSTALLS ######

## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite()

#Required packages
#biocLite(c("GenomicRanges", "Biostrings", "Rsamtools", "DECIPHER"))


#### BIOCONDUCTOR INSTALL END #########

#REQUIRED PACKAGES
library(ape)
library(seqinr)
library(stringr)
library(data.table)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
#library(ips)
library(gdata)
library(DECIPHER)

##############################################################################################
#####################  0. Functions and stuff                    #############################
#####################                                            #############################
##############################################################################################

options(stringsAsFactors = FALSE)
options(warn=2) #for debugging warnings in loops

##############################################################################################
#####################  1.Match loci to contigs                   #############################
#####################                                            #############################
##############################################################################################

#This script does the following:
#1. Matches the loci to the contigs, saves them to a new file
#2. Also finds the potential paralogs, removes them, and saves them to a separate file

#Set up directories
threads<-"6" #threads, keep quotes
probe.type<-"loci" #choose loci if the provided file is loci rather than probes
#probe.type<-"probe"
in.dir<-"Processed_Samples"
out.dir<-"assembled_contigs"
contig.save<-"anolis_uce"

#CLUSTER directories
#work.dir<-"/home/c111h652/scratch/Hyloidea_Analysis"
#contig.dir<-"/home/c111h652/scratch/Hyloidea_Analysis/Assembled_Contigs/dipcontigs"
#probe.file<-"/home/c111h652/scratch/Hyloidea_Analysis/Full_Hyloidea_Loci_May23.fa"


#### HOME COMPUTER

#Ranoidea all
#work.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Phylogeny"
#contig.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Phylogeny/Assembled_Contigs/Ranoidea_Seqcap"
#probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Full_Ranoidea_Loci_May23.fa"

#Hyloidea all
work.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Phylogeny"
contig.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Phylogeny/Assembled_Contigs/Hyloidea_Seqcap"
probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Full_Hyloidea_Loci_Jun27.fa"

#Ranoidea
#work.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Seqcap/Post_Processing"
#contig.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Seqcap/Post_Processing/Assembled_Contigs_Final/Ranoidea_reduced"
#probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Reduced_Ranoidea_Loci_May23.fa"

#Hyloidea
#work.dir<-"/Volumes/Armored/Hyloidea_Analysis"
#contig.dir<-"/Volumes/Armored/Hyloidea_Analysis/Assembled_Contigs"
#probe.file<-"/Volumes/Armored/Hyloidea_Analysis/Final_Hyloidea_Baits_Jun24.fa"
#probe.file<-"/Volumes/Armored/Hyloidea_Analysis/Full_Hyloidea_Loci_Jun27.fa"

#work.dir<-"/Volumes/Armored/Anura_Seqcap"
#contig.dir<-"/Volumes/Armored/Anura_Seqcap/Assembled_Contigs_Final/Hyloidea_broad_24"
#probe.file<-"/Volumes/Armored/Hyloidea_Analysis/Full_Hyloidea_Loci_Jun27.fa"

#Microhylidae
#work.dir<-"/Users/chutter/Dropbox/Research/WIP/Microhylidae_SeqCap"
#contig.dir<-"/Users/chutter/Dropbox/Research/WIP/Microhylidae_SeqCap/Assembled_Contigs_Final"
#probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Full_Ranoidea_Loci_May23.fa"

#Anolis UCE
#work.dir<-"/Users/chutter/Dropbox/Research/WIP/Anolis_Biogeography"
#contig.dir<-"/Users/chutter/Dropbox/Research/WIP/Anolis_Biogeography/Assembled_Contigs"
#probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/uce-5k-probes.fasta"
#probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Hutter_uce5k_loci.fa"

#load in probe file
probe.loci<-scanFa(FaFile(probe.file))

setwd(work.dir)
if (file.exists(in.dir) == F){ dir.create(in.dir) }

#Make blast database for the probe loci
system(paste("makeblastdb -in ", probe.file, " -parse_seqids -dbtype nucl ",
             " -out ", work.dir, "/blast_db", sep = ""))

#headers for the blast db
headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
           "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

#gets lists of directories and files with sample names
setwd(contig.dir)
file.names<-list.files(pattern = "", full.names = F, recursive = T)

########################################################################  
# STEP 1
# Matches contigs to target loci
########################################################################

#Matching and processing for each sample
for (i in 1:length(file.names)){
  
  #Sets up working directories for each species
  sample<-gsub(pattern = ".fa$", replacement = "", x = file.names[i])
  species.dir<-paste(work.dir, "/", in.dir, "/", sample, sep = "")
  #Creates species directory if none exists
  if (file.exists(species.dir) == F){ dir.create(species.dir) }
  if (file.exists(paste(species.dir, "/", out.dir, sep = "")) == F) {dir.create(paste(species.dir, "/", out.dir, sep = "")) }
  setwd(paste(species.dir, "/", out.dir, sep = ""))
  
  # DEDUPE almost exact duplicate removal
  system(paste("dedupe.sh in=",contig.dir, "/", file.names[i], " ordered=t overwrite=true ",
               " out=", sample, "_dd.fa", " minidentity=97", sep = ""), ignore.stderr = T)
  
  #Matches samples to loci
  system(paste("blastn -task dc-megablast -db ", work.dir, "/blast_db",
               " -query ", sample, "_dd.fa -out ", sample, "_match.txt", 
               " -outfmt 6 -num_threads ", threads, sep = ""))
  
  # system(paste("mpirun pblat -threads=", threads, " ", probe.file, 
  #              " ", sample, "_dd.fa -tileSize=8 -minIdentity=60", 
  #              " -noHead -out=pslx ", sample, ".pslx", sep = "")) 
  
  #Need to load in transcriptome for each species and take the matching transcripts to the database
  #match.data<-fread(paste(sample, ".pslx", sep =""), sep = "\t", header = F, stringsAsFactors = FALSE)
  match.data<-fread(paste(sample, "_match.txt", sep =""), sep = "\t", header = F, stringsAsFactors = FALSE)
  setnames(match.data, headers)
  
  if (probe.type == "probe"){
    match.data$tName<-gsub(pattern = "_p.*", replacement = "", x = match.data$tName) 
    match.data$tName<-gsub(pattern = "_.*", replacement = "", x = match.data$tName)
    names(probe.loci)<-gsub(pattern = "_.*", replacement = "", x = names(probe.loci))
  }#end probe if
  
  ##############################################################
  # Part A: Fixes up the match database
  ##############################################################
  
  #Gets rid of very poor matches
  filt.data<-match.data[match.data$matches > 40,]
  filt.data<-filt.data[filt.data$evalue <= 0.05,]
  
  #Fixes direction and adds into data
  filt.data[,qDir:= as.character("0")]
  #Finds out if they are overlapping
  for (k in 1:nrow(filt.data)){
    if (filt.data$tStart[k] > filt.data$tEnd[k]){
      filt.data$qDir[k]<-"-"
      new.start<-min(filt.data$tStart[k], filt.data$tEnd[k])
      new.end<-max(filt.data$tStart[k], filt.data$tEnd[k])
      filt.data$tStart[k]<-new.start
      filt.data$tEnd[k]<-new.end
    } else { filt.data$qDir[k]<-"+" }
  }#end k loop
  
  #Get the sizes from the contig name
  contigs<-scanFa(FaFile(paste(sample, "_dd.fa$", sep = "")))
  new.qsize<-width(contigs)[pmatch(filt.data$qName, names(contigs), duplicates.ok = T)]
  filt.data[,qSize:=as.numeric(new.qsize)]
  #Gets the sizes from the probes
  new.tsize<-width(probe.loci)[pmatch(filt.data$tName, names(probe.loci), duplicates.ok = T)]
  filt.data[,tSize:=as.numeric(new.tsize)]
  #Removes matches that barely match to the full contig
  filt.data<-filt.data[filt.data$matches >= filt.data$tSize *0.20,]
  
  #Fixes duplicate same contigs matching to the same locus (2 entries)
  contig.names<-unique(filt.data[duplicated(filt.data$qName) == T,]$qName)
  save.match<-c()
  for (j in 1:length(contig.names)) {
    sub.match<-filt.data[filt.data$qName %in% contig.names[j],]
    #Skips if 1 row
    if (nrow(sub.match) == 1){ next }

    temp.loci<-unique(sub.match$tName)
    for (k in 1:length(temp.loci)){
      temp.match<-sub.match[sub.match$tName %in% temp.loci[k],]
      #Skips if only 1
      if (nrow(temp.match) == 1){ next }
      #Saves to remove if its a poorly assembled repeat area
      
      temp.match$qStart[1]<-min(temp.match$qStart)
      temp.match$qEnd[1]<-max(temp.match$qEnd)
      temp.match$tStart[1]<-min(temp.match$tStart)
      temp.match$tEnd[1]<-max(temp.match$tEnd)
      temp.match$bitscore[1]<-max(temp.match$bitscore)
      temp.match$evalue[1]<-mean(temp.match$evalue)
      temp.match$pident[1]<-mean(temp.match$pident)
      temp.match$matches[1]<-temp.match$qEnd[1]-temp.match$qStart[1]
      save.match<-rbind(save.match, temp.match[1,])
    }#end k loop
  }#end j loop
  
  new.data<-filt.data[!filt.data$qName %in% save.match$qName,]
  new.data<-rbind(new.data, save.match)

  ##############################################################
  # Part B: Finds different contigs that match to the same probe locus
  ##############################################################
  
  #next go through all the remaining duplicates
  loci.names<-unique(new.data[duplicated(new.data$tName) == T,]$tName)
  fix.seq<-DNAStringSet()
  save.paralog<-c()
  del.loci<-c()
  for (j in 1:length(loci.names)){
    #subsets data
    sub.match<-new.data[new.data$tName %in% loci.names[j],]
    sub.match<-sub.match[order(sub.match$tStart)]
    if (length(unique(sub.match$qName)) == 1 && length(unique(sub.match$tName)) == 1){ next }
    
    #removes bad things
    contig.per<-sub.match$matches/sub.match$tSize
    sub.match<-sub.match[contig.per >= 0.15,]
    if (nrow(sub.match) == 0){ 
      save.paralog<-append(save.paralog, loci.names[j])
      next }
    
    #Keep if they match to same contig, then not a paralog
    if (length(unique(sub.match$qName)) == 1){
      save.contig<-contigs[names(contigs) %in% sub.match$qName[1]]
      names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
      fix.seq<-append(fix.seq, save.contig)
      next
    }
    
    #If the number is negative then problem!
    hit.para<-0
    for (k in 1:(nrow(sub.match)-1)){
      if (sub.match$tStart[k+1]-sub.match$tEnd[k] < -20){ hit.para<-1 }
    }  
    
    #If there are overlaps
    if (hit.para == 1){
      save.paralog<-append(save.paralog, sub.match$tName[1])
      next
    }#end if
    
    #If there is a single match above 50% keep it its good
    per.match<-sub.match[(sub.match$matches/sub.match$tSize) >= 0.9,]
    if (nrow(per.match) == 1){ 
      save.contig<-contigs[names(contigs) %in% sub.match$qName[1]]
      names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
      fix.seq<-append(fix.seq, save.contig)
      next
    }
    
    #If the number is negative then problem!
    row.rem<-c()
    for (k in 1:(nrow(sub.match)-1)){
      if (sub.match$tStart[k+1]-sub.match$tEnd[k] < -(sub.match$matches[k]*0.20) ){ 
        row.rem<-append(row.rem, sub.match$qName[k])
        row.rem<-append(row.rem, sub.match$qName[k+1])
      }#end if
    }  #end k
    
    #Removes overlapping matches
    red.match<-sub.match[!sub.match$qName %in% row.rem,]
    if (nrow(red.match) == 0){  
      #If there is a single match above 50% keep it its good
      del.loci<-append(del.loci, row.rem)
      per.match<-sub.match[(sub.match$matches/sub.match$tSize) >= 0.5,]
      if (nrow(per.match) == 1){ 
        save.contig<-contigs[names(contigs) %in% per.match$qName[1]]
        names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
        fix.seq<-append(fix.seq, save.contig)
        next
      }
      next
    }#end if statement
    
    #Looks for a pretty high bit score
    b.match<-sub.match[sub.match$bitscore/max(sub.match$bitscore) >= 0.3,]
    if (nrow(b.match) == 1){ 
      save.contig<-contigs[names(contigs) %in% b.match$qName[1]]
      names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
      fix.seq<-append(fix.seq, save.contig)
      next
    }
    
    #Checks if the two pieces may overlap
    hit.over<-0
    overlap.count<-0
    for (k in 1:(nrow(sub.match)-1)){
      overlap.end<-sub.match$tEnd[k]+(sub.match$qSize[k]-sub.match$qEnd[k])
      overlap.start<-sub.match$tStart[k+1]-sub.match$qStart[k+1]
      
      #checks if they overlap
      if (overlap.start < overlap.end) {
        hit.over<-1
        break
      }
    }#end k loop
      
    #If they overlap, then cut out edge intron as this is probably two exons    
    if (hit.over == 1){
      #Cuts the node apart and saves separately
      sub.match$qStart[1]<-as.numeric(1)
      sub.match$qEnd[nrow(sub.match)]<-sub.match$qSize[nrow(sub.match)]
      #If it ends up with a negative start
  
      #Collects new sequence fragments
      spp.seq<-contigs[names(contigs) %in% sub.match$qName]
      spp.seq<-spp.seq[pmatch(sub.match$qName, names(spp.seq))]
      
      new.seq<-DNAStringSet()
      for (k in 1:length(spp.seq)){ new.seq<-append(new.seq, subseq(x = spp.seq[k], start = sub.match$qStart[k], end = sub.match$qEnd[k]) ) }  
      
      #Combine new sequence
      save.contig<-DNAStringSet(paste(as.character(new.seq), collapse = "", sep = "") )
      names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
      fix.seq<-append(fix.seq, save.contig)
      next
    }#end if 
    
    #Sets up the new contig location
    #Cuts the node apart and saves separately
    sub.match$tEnd<-sub.match$tEnd+(sub.match$qSize-sub.match$qEnd)
    sub.contigs<-contigs[names(contigs) %in% sub.match$qName]
    
    join.contigs<-DNAStringSet()
    for (k in 1:(nrow(sub.match)-1)){
      join.contigs<-append(join.contigs, sub.contigs[k])
      n.pad<-sub.match$tStart[k+1]-sub.match$tEnd[k]
      join.contigs<-append(join.contigs, DNAStringSet(paste(rep("N", n.pad), collapse = "", sep = "")) )
    }
    join.contigs<-append(join.contigs, sub.contigs[length(sub.contigs)])
    save.contig<-DNAStringSet(paste(as.character(join.contigs), collapse = "", sep = "") )
    
    #Saves final sequence
    names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
    fix.seq<-append(fix.seq, save.contig)
  }# end j loop
  
  ##############################################################
  # Part C: Finds the same contig that match to different probe locus
  ##############################################################
  #gets rid of double dups
  temp.data<-new.data[!new.data$tName %in% save.paralog,]
  temp.data<-temp.data[!temp.data$qName %in% del.loci,]
  dd.data<-temp.data[!temp.data$tName %in% gsub("_\\|_.*", "", names(fix.seq)),]
  #Find contigs that match to multiple loci
  temp.q<-dd.data[duplicated(dd.data$qName) == T,]$qName
  #Pulls out the double match data
  dup.match<-dd.data[dd.data$qName %in% temp.q,]
  #Gets names of things that are surely good because the contig is in the matches twice
  dup.data<-dup.match[order(dup.match$qName)]
  
  #Loops through each potential duplicate
  dup.loci<-unique(dup.data$qName)
  for (j in 1:length(dup.loci)){
    #pulls out data that matches to multiple contigs
    sub.data<-dup.data[dup.data$qName %in% dup.loci[j],]
    sub.data<-sub.data[order(sub.data$qStart)]

    #Saves them if it is split up across the same locus
    if (length(unique(sub.data$tName)) == 1 && length(unique(sub.data$qName)) == 1){
      spp.seq<-contigs[names(contigs) %in% sub.data$qName]
      names(spp.seq)<-paste(sub.data$tName[1], "_|_", sample, sep ="")
      fix.seq<-append(fix.seq, spp.seq)
      next
    }
    
    #Cuts the node apart and saves separately
    sub.data$qStart<-sub.data$qStart-sub.data$tStart
    #If it ends up with a negative start
    sub.data$qStart[sub.data$qStart <= 0]<-1
    #Fixes ends
    sub.data$qEnd<-sub.data$qEnd+(sub.data$tSize-sub.data$tEnd)
    
    #Fixes if the contig is smaller than the full target locus
    sub.data$qEnd[sub.data$qEnd >= sub.data$qSize]<-sub.data$qSize[1]
    
    starts<-c()
    ends<-c()
    starts[1]<-1
    for (k in 1:(nrow(sub.data)-1)){
      ends[k]<-sub.data$qEnd[k]+floor((sub.data$qStart[k+1]-sub.data$qEnd[k])/2)
      starts[k+1]<-ends[k]+1
    } #end k loop
    ends<-append(ends, sub.data$qSize[1])
    
    #Looks for overlapping contigs
    tmp<-ends-starts
    if(length(tmp[tmp < 0 ]) != 0){ 
      save.paralog<-append(save.paralog, sub.match$tName)
      next
    }
    #Collects new sequence fragments
    spp.seq<-contigs[names(contigs) %in% sub.data$qName]
    new.seq<-DNAStringSet()
    for (k in 1:length(starts)){ new.seq<-append(new.seq, subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }  
    
    #renames and saves
    names(new.seq)<-paste(sub.data$tName, "_|_", sample, sep ="")
    fix.seq<-append(fix.seq, new.seq)
  } #end j loop
  
  #Writes the base loci
  temp.data<-dd.data[!dd.data$tName %in% gsub("_\\|_.*", "", names(fix.seq)),]
  base.data<-temp.data[!temp.data$tName %in% save.paralog,]
  base.loci<-contigs[names(contigs) %in% base.data$qName]
  sort.data<-base.data[match(names(base.loci), base.data$qName),]
  #Name and finalize
  names(base.loci)<-paste(sort.data$tName, "_|_", sample, sep ="")
  fin.loci<-append(base.loci, fix.seq)
  fin.loci<-fin.loci[width(fin.loci) >= 60]
  
  #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
  temp<-fin.loci[duplicated(names(fin.loci)) == T]
  if(length(temp) != 0){ stop("DUPLICATE FOUND") }
  
  #Finds probes that match to two or more contigs
  final.loci<-as.list(as.character(fin.loci))
  write.fasta(sequences = final.loci, names = names(final.loci), 
              paste(sample, "_orthologs.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  paralog.names<-new.data[new.data$tName %in% save.paralog,]
  paralog.contigs<-contigs[names(contigs) %in% paralog.names$qName]
  write.loci<-as.list(as.character(paralog.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(sample, "_paralogs.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  system(paste("rm ", sample, "_dd.fa ", sep = ""))  
  system(paste("rm ", sample, "_match.txt", sep = ""))  
  print(paste(sample, " probe matching complete. ", length(final.loci), " found!"))
  
}# end i loop

########################################################################  
# STEP 2
# Output a file of summary stats
########################################################################

#Sets up data summary
setwd(work.dir)
system(paste("rm blast_db*", sep = ""))  

header.data<-c("Sample", "noContigs", "orthoLoci", "paraLoci", "N50", "N90", "minLen", "maxLen", "meanLen", "medianLen")   
samples<-gsub(".fa$", "", file.names)
prelim.data<-data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
setnames(prelim.data, header.data)
prelim.data[, Sample:=as.character(samples)]
merge.contigs<-DNAStringSet()

#Cycles through each assembly run and assesses each
for (i in 1:length(samples)){
  
  #Gets raw contigs to count them
  setwd(paste(work.dir, "/", in.dir, "/", samples[i], "/", out.dir, sep = ""))
 
  #Gets length of raw contigs
  contigs<-scanFa(FaFile(paste(contig.dir, "/", samples[i], ".fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("noContigs", header.data), value = length(contigs) )
  
  #Gets length of raw contigs
  og<-scanFa(FaFile(paste(samples[i], "_orthologs.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("orthoLoci", header.data), value = length(og) )
  
  #Gets length of raw contigs
  para<-scanFa(FaFile(paste(samples[i], "_paralogs.fa", sep = "")))
  temp.para<-unique(names(para))
  set(prelim.data, i = match(samples[i], samples), j = match("paraLoci", header.data), value = length(temp.para) )
  
  #Summarizes basic data
  len.sorted <- rev(sort(width(og)))
  N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1]
  N90 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.9][1]
  
  set(prelim.data, i =  match(samples[i], samples), j = match("N50", header.data), value = N50 )
  set(prelim.data, i =  match(samples[i], samples), j = match("N90", header.data), value = N90 )
  set(prelim.data, i =  match(samples[i], samples), j = match("minLen", header.data), value = min(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("maxLen", header.data), value = max(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("meanLen", header.data), value = mean(width(og)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("medianLen", header.data), value = median(width(og)) )
  
  merge.contigs<-append(merge.contigs, og)
  
}#End loop for things

#Saves combined, final dataset
setwd(work.dir)
write.csv(prelim.data, file = paste(contig.save, "-sample-assessment.csv", sep = ""))

final.loci<-as.list(as.character(merge.contigs))
write.fasta(sequences = final.loci, names = names(final.loci), paste(contig.save, "_contigs.fa", sep = ""), nbchar = 1000000, as.string = T)

#### END SCRIPT 
# New files will be located in each species folder within 'Processed_Samples', in the folder assembled_loci
