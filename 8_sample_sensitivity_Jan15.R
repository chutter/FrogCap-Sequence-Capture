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

options(stringsAsFactors = FALSE)
options(warn=2) #for debugging warnings in loops


###############################################################################
###############################################################################
######################Step 1 pull out loci and align  #########################
###############################################################################
###############################################################################

work.dir<-"/Volumes/Armored/Test_dataset"
contig.dir<-"/Volumes/Armored/Test_dataset/dipcontigs"
out.dir<-"Loci_Matching"
threads<-"6" #threads, keep quotes

#First is for UCE, second is for Ranoidea
probe.file<-"/Volumes/Armored/Ranoidea_Analysis/Final_Datasets_2017/Final_Ranoidea_Loci_Jun3.fa"

#gets lists of directories and files with sample names
setwd(contig.dir)
file.names<-list.files(pattern = "", full.names = F, recursive = T)

#Look for paralogs and save
dir.create(paste(work.dir, "/", out.dir, sep = ""))
setwd(paste(work.dir, "/", out.dir, sep = ""))

#PSLX headers
headers<-c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
           "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSize", "qStarts", "tStarts", "qSeq", "tSeq")

########################################################################  
# STEP 1
# Matches contigs to target loci
########################################################################

#Matching and processing for each sample
for (i in 2:length(file.names)){
  #Gets the sample name
  samples<-gsub(pattern = ".fa$", replacement = "", x = file.names[i])
  
  #creates sample directory
  samp.dir<-paste(work.dir, "/", out.dir, "/", samples, sep = "")
  dir.create(samp.dir)
  setwd(samp.dir)
  
  # DEDUPE almost exact duplicate removal
  system(paste("dedupe.sh in=",contig.dir, "/", file.names[i], " ordered=t ",
               "out=", samples, "_dd.fa", " minidentity=97", sep = ""))
  
  #Matches samples to loci
  system(paste("mpirun pblat -threads=", threads, " ", probe.file, 
               " ", samp.dir, "/", samples, "_dd.fa -tileSize=8 -minIdentity=60", 
               " -noHead -out=pslx ", samples, ".pslx", sep = "")) 
  
  #Need to load in transcriptome for each species and take the matching transcripts to the database
  match.data<-fread(paste(samples, ".pslx", sep =""), sep = "\t", header = F, stringsAsFactors = FALSE)
  setnames(match.data, headers)
  
  #Remove poor, small matches
  match.data<-match.data[match.data$matches > 60,]
  filt.data<-match.data[(match.data$tEnd-match.data$tStart) >= (match.data$tSize * 0.30),]
  
  #Step 2, fixes this problem
  #Contig = Locus (both same)
  #Contig = Locus (both same)
  
  #Find contigs that match to multiple loci
  temp.t<-filt.data[duplicated(filt.data$tName) == T,]$tName
  #Pulls out the double match data
  t.match<-filt.data[filt.data$tName %in% temp.t,]
  #Gets names of things that are surely good because the contig is in the matches twice
  t.match<-t.match[order(t.match$tName)]
  
  #gets locus names
  loci.names<-unique(t.match$tName)
  to.keep<-c()
  #loops through each locus and assesses if they are matching to same contig
  for (j in 1:length(loci.names)){
    #subsets data
    sub.match<-t.match[t.match$tName %in% loci.names[j],]
    
    #Removes if they match to same contig, then not a paralog
    temp.dup<-unique(sub.match$qName)
    if (length(temp.dup) == 1){
      to.keep<-append(to.keep, loci.names[j])
      next
    } #end if
  }# end j loop
  
  #gets rid of double dups
  t.data<-t.match[t.match$tName %in% to.keep,]
  temp.data<-filt.data[!filt.data$tName %in% to.keep,]
  t.data<-t.data[duplicated(t.data$tName) != T,]
  dd.data<-rbind(temp.data, t.data)
  
  #Find contigs that match to multiple loci
  temp.q<-dd.data[duplicated(dd.data$qName) == T,]$qName
  #Pulls out the double match data
  dup.match<-dd.data[dd.data$qName %in% temp.q,]
  #Gets names of things that are surely good because the contig is in the matches twice
  dup.data<-dup.match[order(dup.match$qName)]
  
  #Loops through each locus and fixes
  fa<-FaFile(paste(samples, "_dd.fa", sep = "")) 
  contigs<-scanFa(fa)
  names(contigs)<-gsub(pattern = " .*", replacement = "", x = names(contigs))
  
  dup.loci<-unique(dup.data$qName)
  fixed.loci<-DNAStringSet()
  for (j in 1:length(dup.loci)){
    #pulls out data that matches to multiple contigs
    sub.data<-dup.data[dup.data$qName %in% dup.loci[j],]
    sub.data<-sub.data[order(sub.data$qStart)]
    
    #Saves them if it is split up across the same locus
    if (length(unique(sub.data$tName)) == 1){
      spp.seq<-contigs[names(contigs) %in% sub.data$qName]
      names(spp.seq)<-paste(sub.data$tName[1], "_|_", samples, sep ="")
      fixed.loci<-append(fixed.loci, spp.seq)
      next
    }
    
    #Removes duplicates  
    sub.data<-sub.data[duplicated(sub.data$tName) !=T,]
    
    #Remove too long things
    sub.data<-sub.data[(sum(sub.data$matches)/sub.data$tSize[1]) > 0.30, ]
    sub.data<-sub.data[sub.data$qBaseInsert <= 100, ]
    if (nrow(sub.data) == 0) { next }
    
    #Saves them if it is split up across the same locus
    if (nrow(sub.data) == 1){
      spp.seq<-contigs[names(contigs) %in% sub.data$qName]
      names(spp.seq)<-paste(sub.data$tName[1], "_|_", samples, sep ="")
      fixed.loci<-append(fixed.loci, spp.seq)
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
    
    tmp<-ends-starts
    if(length(tmp[tmp < 0 ]) != 0){ next }
    #Collects new sequence fragments
    spp.seq<-contigs[names(contigs) %in% sub.data$qName]
    new.seq<-DNAStringSet()
    for (k in 1:length(starts)){ new.seq<-append(new.seq, subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }  
    
    #renames and saves
    names(new.seq)<-paste(sub.data$tName, "_|_", samples, sep ="")
    fixed.loci<-append(fixed.loci, new.seq)

  } #end j loop
  
  #Writes the base loci
  base.data<-dd.data[!dd.data$qName %in% dup.loci, ]
  base.loci<-contigs[names(contigs) %in% base.data$qName]
  sort.data<-base.data[match(names(base.loci), base.data$qName),]
  names(base.loci)<-paste(sort.data$tName, "_|_", samples, sep ="")
  fin.loci<-append(base.loci, fixed.loci)
  
  #Finds probes that match to two or more contigs
  temp.p<-fin.loci[duplicated(names(fin.loci)) == T]
  #Pulls out the double match data
  par.loci<-fin.loci[names(fin.loci) %in% names(temp.p)]
  
  #Writes the paralog loci
  final.loci<-as.list(as.character(par.loci))
  write.fasta(sequences = final.loci, names = names(final.loci), 
              paste(samples, "_paralog.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  red.loci<-fin.loci[!names(fin.loci) %in% names(par.loci)]  
  final.loci<-as.list(as.character(red.loci))
  write.fasta(sequences = final.loci, names = names(final.loci), 
              paste(samples, "_base.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  system(paste("rm ", samples, "_dd.fa ", samples, ".pslx", sep = ""))  
  print(paste(samples, " probe matching complete"))
  
}# end i loop

#Goes to working directory and pulls out sample names
setwd(paste(work.dir, "/", out.dir, sep = ""))
sample.names<-list.dirs(path = ".", full.names = F, recursive = F)

#Prelim data
header.data<-c("Sample", "Contigs", "lociMatches", "noParalog", "perMatch", "N50", "N90", 
               "minLen", "maxLen", "meanLen", "medianLen")
prelim.data<-data.table(matrix(as.double(0), nrow = length(sample.names), ncol = length(header.data)))
setnames(prelim.data, header.data)
prelim.data[, Sample:=as.character(sample.names)]

#Cycles through each assembly run and assesses each
for (i in 1:length(sample.names)){
  
  #creates sample directory
  samp.dir<-paste(work.dir, "/", out.dir, "/", sample.names[i], sep = "")
  setwd(samp.dir)
  
  #Matches samples to loci
  system(paste("mpirun pblat -threads=", threads, " ", probe.file, 
               " ", samp.dir, "/", sample.names[i], "_base.fa -tileSize=8 -minIdentity=60", 
               " -noHead -out=pslx ", sample.names[i], ".pslx", sep = "")) 
  
  #Need to load in transcriptome for each species and take the matching transcripts to the database
  match.data<-fread(paste(sample.names[i], ".pslx", sep =""), sep = "\t", header = F, stringsAsFactors = FALSE)
  setnames(match.data, headers)
  
  #Gets names of things that are surely good because the contig is in the matches twice
  match.data<-match.data[order(match.data$tName)]
  
  #gets locus names
  loci.names<-unique(match.data$tName)
  #loops through each locus and assesses if they are matching to same contig
  spp.spec<-c()
  for (j in 1:length(loci.names)){
    #subsets data
    sub.match<-match.data[match.data$tName %in% loci.names[j],]
    spp.spec<-append(spp.spec, (sub.match$matches+sub.match$misMatches)/sub.match$tSize)
  }# end j loop

  print(paste(sample.names[i], mean(spp.spec)))
  
} #end loop
  
  
  
  #Gets length of raw contigs
  fa<-FaFile(paste(sample.names[i], ".fa", sep = "")) 
  set(prelim.data, i = match(sample.names[i], sample.names), j = match("Contigs", header.data), value = length(names(scanFa(fa))) )
  
  #Gets assembly metrics for matching contigs
  samp.dir<-paste(work.dir, "/", out.dir, "/", sample.names[i], sep = "")
  setwd(samp.dir)
  
  #Reads in paralog
  par<-FaFile(paste(sample.names[i], "_paralog.fa", sep = "")) 
  
  #Loops through each locus and fixes
  loc<-FaFile(paste(sample.names[i], "_base.fa", sep = "")) 
  contigs<-scanFa(loc)
  
  len.sorted <- rev(sort(width(contigs)))
  N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1]
  N90 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.9][1]
  
  #Summarizes basic data
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("lociMatches", header.data), value = length(contigs) )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("noParalog", header.data), value = length(names(scanFa(par))) )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("perMatch", header.data), value =  length(contigs)/5000 )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("N50", header.data), value = N50 )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("N90", header.data), value = N90 )
  
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("minLen", header.data), value = min(width(contigs)) )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("maxLen", header.data), value = max(width(contigs)) )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("meanLen", header.data), value = mean(width(contigs)) )
  set(prelim.data, i =  match(sample.names[i], sample.names), j = match("medianLen", header.data), value = median(width(contigs)) )
  
}#End loop for things

