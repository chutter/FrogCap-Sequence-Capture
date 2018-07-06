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
library(gdata)
library(ShortRead)

options(stringsAsFactors = FALSE)
options(warn=2) #for debugging warnings in loops

spades<-function(spades.path, read.pe1, read.pe2, read.pes){
  
 # read.pe1<-"a_read1.fastq"
 # read.pe2<-"a_read2.fastq"
 # read.pes<-"a_read3.fastq"
  print("Spades run begin.")
  k<-c(9,13,21,33,55,77,99,127)
  k.val<-paste(k, collapse = ",")
  system(paste(spades.path, " --pe1-1 ", read.pe1, " --pe1-2 ", read.pe2, " --pe1-s ", read.pes,
               " -o spades -k ",k.val," --careful --expect-gaps --hap-assembly ",
               " -t ", threads, " -m ", mem, sep = ""), ignore.stdout = T) 
  
  #Does quality checks in case spades fails for some reason
  if (file.exists("spades/spades/contigs.fasta") == T){
    if (file.info("spades/spades/contigs.fasta")$size == 0){ system("rm spades/spades/contigs.fasta")}
  }  
  phred.off<-""
  while (file.exists("spades/spades/contigs.fasta") == F){
    #First checks to see if the error is due to some phred score thing
    temp.dirs<-list.dirs("spades/spades", full.names = F, recursive = F)
    k.pres<-grep("K", temp.dirs)
    #First checks for weird phred error
    if (length(k.pres) == 0){
      phred.off<-"--phred-offset 33"
      system(paste(spades.path, " --pe1-1 ", read.pe1, " --pe1-2 ", read.pe2, " --pe1-s ", read.pes,
                   " -o spades -k ",k.val," --careful --expect-gaps --hap-assembly ", phred.off, 
                   " -t ", threads, " -m ", mem, sep = ""), ignore.stdout = T) 
    } #end if
    
    #Next moves on to see if missing k folder caused error
    temp.dirs<-list.dirs("spades/spades", full.names = F, recursive = F)
    k.pres<-temp.dirs[grep("K", temp.dirs)]
    
    #subtract Ks until it works, exits with no data if it does not
    k<-k[-(length(k.pres))]
    if (length(k) == 0) { 
      print(paste("k-mer values all used up, cannot assemble!"))
      if (file.exists("spades") == T) { system("rm -r spades") }
      return(NULL)
    }
    
    #Tries the newly subtracted k values
    k.val<-paste(k, collapse = ",")
    system("rm -r spades")
    system(paste(spades.path, " --pe1-1 ", read.pe1, " --pe1-2 ", read.pe2, " --pe1-s ", read.pes,
                 " -o spades -k ",k.val," --careful --expect-gaps --hap-assembly ", phred.off, 
                 " -t ", threads, " -m ", mem, sep = ""), ignore.stdout = T) 
    
    if (file.exists("spades/spades/contigs.fasta") == T){
      if (file.info("spades/spades/contigs.fasta")$size == 0){ system("rm spades/spades/contigs.fasta")}
    }
  }#end while
  
  #If the k-mers are all run out, therefore nothing can be assembled
  if (length(k) == 0) { 
    print(paste("k-mer values all used up, cannot assemble!"))
    if (file.exists("spades") == T) { system("rm -r spades") }
    return(NULL)
  }# end if 
  
  if (file.exists("spades/dipspades/consensus_contigs.fasta") != T){
    best.contig<-scanFa(FaFile("spades/spades/contigs.fasta"))   # loads up fasta file
  } else {
    sp.contig<-scanFa(FaFile("spades/spades/contigs.fasta"))   # loads up fasta file
    ds.contig<-scanFa(FaFile("spades/dipspades/consensus_contigs.fasta"))   # loads up fasta file
    cn.contig<-scanFa(FaFile("spades/dipspades/possibly_conservative_regions.fasta"))  
  
    if (max(width(cn.contig)) > max(width(ds.contig))) { best.contig<-cn.contig 
    } else if (max(width(ds.contig)) >= max(width(sp.contig))) { best.contig<-ds.contig 
    } else { best.contig<-sp.contig }
  }#end else
    system("rm -r spades")
    return(best.contig)
}

cap3<-function(input.contigs, min.contig.length){
  
 # input.contigs<-contigs
  #min.contig.length<-50
  #Writes contigs for cap3
  write.loci<-as.list(as.character(input.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              "temp_seed.fa", nbchar = 1000000, as.string = T)
  
  #runs cap3 to merge similar contigs (pull only clustered contigs out?)
  system(paste("cap3 temp_seed.fa -z 1 -o 16 -e 11 -s 251", " > ", 
               "log.fa.cap.txt", sep = "")) 
  
  #Reads in results files
  temp.assembled<-scanFa(FaFile(paste("temp_seed.fa.cap.contigs", sep = "")))
  temp.singlets<-scanFa(FaFile(paste("temp_seed.fa.cap.singlets", sep = "")))
  keep.singlets<-temp.singlets[width(temp.singlets) >= min.contig.length]
  if (length(temp.assembled) != 0) { names(temp.assembled)<-paste("cap3_assembled_", seq(1:length(temp.assembled)), sep = "") }
  if (length(keep.singlets) != 0) {names(keep.singlets)<-paste("cap3_singlets_", seq(1:length(keep.singlets)), sep = "") }
  cap3.output<-append(temp.assembled, keep.singlets)
  
  #Get cap3 files and deletes
  cap.files<-list.files(pattern = "", full.names = F, recursive = F)
  cap.remove<-cap.files[grep(pattern = paste("fa.cap*.", sep =""), x = cap.files)]
  system(paste("rm ", paste(cap.remove, collapse = " ") ))
  system("rm temp_seed.fa")
  
  return(cap3.output)
  
}

mafft<-function(unaligned.contigs){
  
  unaligned.contigs<-save.align
  
  #Align contigs with reference to see which matches
  final.loci<-as.list(as.character(unaligned.contigs))
  write.fasta(sequences = final.loci, names = names(final.loci), 
              "mafft_temp.fa", nbchar = 1000000, as.string = T)
  
  #Runs MAFFT to align
  system(paste("mafft --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123",
               " --thread ", threads, " mafft_temp.fa > temp_align.fa", sep = ""))
  
  alignment<-scanFa(FaFile("temp_align.fa", sep = ""))   # loads up fasta file
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  
  
}

pairwise.inf.sites<-function(x, y) {
  #Sets up data, puts ref seq first
  new.align<-x
  new.align[new.align == "n"]<-"-"
  new.align[is.na(new.align) == T]<-"-"
  ref<-new.align[rownames(new.align) == y,]
  summary.data<-c()
  all.pars<-c()
  all.over<-c()
  for (z in 1:nrow(new.align)) {
    #Site counter
    pars<-0
    overlap<-0
    tar<-new.align[z,]
    combined<-matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,]<-ref
    combined[2,]<-tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col<-vector("character", length = nrow(combined))
      seq.col<-combined[,k]
      #not equal to -
      f.char<-seq.col[seq.col != '-'] 
      #don't count missing seq
      if (length(f.char) <= 1) { next }
      
      if (length(f.char) >= 2){
        overlap<-overlap+1
        if (f.char[1] != f.char [2]) { pars<-pars+1 }
      }#end if
    }#ends informative sites loop
    all.pars<-append(all.pars, pars)
    all.over<-append(all.over, overlap)
  }# ends seq loop
  #Summarizes and returns data
  summary.data<-all.pars/all.over
  summary.data[is.nan(summary.data)]<-0
  names(summary.data)<-rownames(new.align)
  return(summary.data)
}

trim.ends<-function (x, min.n.seq = 4){
  if (!inherits(x, "DNAbin")) {
    stop("'x' is not of class 'DNAbin'")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a matrix")
  }
  replaceWithN <- function(x) {
    id <- x == as.raw(4)
    if (length(id) > 0 & any(id[c(1, length(id))])) {
      id <- which(id)
      getIndex <- function(x) {
        for (i in seq_along(id) - 1) {
          if (any(id[1:(i + 1)] != (1:(i + 1)))) 
            break
        }
        id <- rev(id)
        jj <- head(id, 1)
        j <- jj - 1
        for (k in seq_along(id)[-1]) {
          if (any(id[1:k] != (jj:j))) 
            break
          j <- j - 1
        }
        j <- j + 1
        id <- c(0:i, j:jj)
        id[id != 0]
      }
      id <- getIndex(id)
      x[id] <- as.raw(240)
    }
    return(x)
  }
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b) {
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if (max(m) < min.n.seq) 
    stop("alignment contains less sequences then required")
  m <- range(which(m >= min.n.seq))
  m <- seq(from = m[1], to = m[2])
  x <- x[, m]
  x
}

write.phy<-function (x, file = "", interleave = FALSE, strict = FALSE){
  str2cha <- function(x) {
    unlist(strsplit(x, ""))
  }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""), 
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i], 
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames))) 
      cat("WARNING: Truncation of taxon names created", 
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ", 
                                                                  diff[i]), collapse = ""), sep = "")
  }
  if (!interleave) 
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1) 
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }
  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm))) 
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1) 
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1]))) 
      sm <- c(sm, "")
    phy <- c(phy, sm)
  }
  if (file == "") {
    cat(phy, sep = "\n")
  }
  else {
    write(phy, file = file)
  }
}


##############################################################################################
#####################  1.Match loci to contigs                   #############################
#####################                                            #############################
##############################################################################################

#This script does the following:
#1. Matches the loci to the contigs, saves them to a new file
#2. Also finds the potential paralogs, removes them, and saves them to a separate file

#Set up directories
threads<-"6" #threads, keep quotes
mem<-"10"
proc.dir<-"Processed_Samples"
raw.dir<-"assembly_reads" #Directory of reads used for assembly. Shouldn't need to modify unless the assembly reads fucked
spades.path<-"dipspades.py"
work.dir<-"/home/c111h652/scratch/Reduced_Probes"
probe.file<-"/home/c111h652/scratch/Reduced_Probes/Reduced_Ranoidea_Loci_Apr21.fa"

probe.file<-"/Volumes/Armored/Reduced_Probes/Reduced_Ranoidea_Loci_Apr21.fa"
spades.path<-"/usr/local/spades/bin/dipspades.py"
work.dir<-"/Volumes/Armored/Reduced_Probes"

#Sets up the reads
setwd(paste(work.dir, "/", proc.dir, sep = ""))
files<-list.files(path = ".", full.names = F, recursive = T)
samples<-list.dirs(path = ".", full.names = F, recursive = F)
reads<-files[grep(pattern = raw.dir, x = files)]
loci.seqs<-scanFa(FaFile(probe.file))
loci.names<-names(loci.seqs)

#PSLX headers
headers<-c("matches", "misMatches", "repMatches", "nCount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", 
           "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSize", "qStarts", "tStarts", "qSeq", "tSeq")

#################################################################
#Step 1: Gather read data and assemble mitochondrial genomes
#################################################################
for (i in 1:length(samples)){
  
  #Change to main directory
  setwd(paste(work.dir, "/", proc.dir, sep = ""))
  ass.dir<-paste(work.dir, "/", proc.dir, "/", samples[i], "/assembled_loci", sep = "")
  if (file.exists(ass.dir) == F) {dir.create(ass.dir) }
  setwd(ass.dir)
  
  #Gets reads together for this sample
  probe.matched<-scanFa(FaFile(paste(samples[i], "_originals.fa", sep = "")))
  sample.reads<-reads[grep(samples[i], reads)]
  read1<-paste(work.dir, "/", proc.dir, "/", sample.reads[grep("READ1", sample.reads)], sep = "")
  read2<-paste(work.dir, "/", proc.dir, "/", sample.reads[grep("READ2", sample.reads)], sep = "")
  read3<-paste(work.dir, "/", proc.dir, "/", sample.reads[grep("singleton", sample.reads)], sep = "")
  
  #Cycles through each locus and tries to assemble ##### CHECK LONGEST 44
  min.id<-"0.7"
  count<-0
  bigger.contigs<-DNAStringSet()
  new.contigs<-DNAStringSet()
  paralog.contigs<-DNAStringSet()
  for (j in 1:length(loci.names)){
    
    #####################################################
    #Step 1. Set up reference locus and first pool of reads
    #####################################################  
    
    #Saves the reference from the locus files
    temp.locus<-loci.seqs[names(loci.seqs) == loci.names[j]]
    temp.assem<-probe.matched[grep(names(temp.locus), names(probe.matched))]
    if (length(temp.assem) != 0 ){
      if (width(temp.assem) >= width(temp.locus)*.8) { next }
      #break
    }
    
    #Writes if its going to try and reassemble
    #comb.ref<-append(temp.locus, temp.assem)
    final.loci<-as.list(as.character(temp.locus))
    write.fasta(sequences = final.loci, names = names(final.loci), 
                paste(ass.dir, "/reference.fa", sep = ""), nbchar = 1000000, as.string = T)
    
    #####################################################
    #Step 2. While loop iterative baiting for assembly
    #####################################################  
    
    #Sets up while loop empty parameters
    system("touch current_seed.fa a_read1.fastq a_read2.fastq a_read3.fastq")
    m.count<-0
    new.len<-0
    counter<-0
    repeat.counter<-0
    seeding<-T
    no.data<-F
    while (seeding == T){
      
      #####################################################
      #While 1. Sets up each round of read data
      #####################################################
      #Copy new reference to do recursively
      counter<-counter+1
      prev.len<-new.len
      
      #Pick out matching reads to the reference for directions and singletons
      system(paste("bbmap.sh -Xmx", mem, "g ref=reference.fa", " in1=", read1, " in2=", read2, " vslow k=12 minid=",min.id, 
                   " outm1=read1.fastq outm2=read2.fastq", sep = ""), ignore.stderr = T)
      system(paste("bbmap.sh -Xmx", mem, "g ref=reference.fa", " in=", read3, " vslow k=12 minid=", min.id,
                   " outm=read3.fastq", sep = ""), ignore.stderr = T)
      
      #Checks if there are any matching reads. If not, moves on
      temp.count<-sum(file.info("read1.fastq")$size, file.info("read2.fastq")$size, file.info("read3.fastq")$size,
                      file.info("a_read1.fastq")$size, file.info("a_read2.fastq")$size, file.info("a_read3.fastq")$size)
        
      if (temp.count == 0) { 
        system("rm -r ref read1.fastq read2.fastq read3.fastq reference.fa")
        seeding<-F
        no.data<-T
        next
        }
      
      #Adds them onto previously found reads
      system(paste("cat a_read1.fastq read1.fastq >> t_read1.fastq"))
      system(paste("cat a_read2.fastq read2.fastq >> t_read2.fastq"))
      system(paste("cat a_read3.fastq read3.fastq >> t_read3.fastq"))
      system("rm read1.fastq read2.fastq read3.fastq a_read1.fastq a_read2.fastq a_read3.fastq")
      system(paste("mv t_read1.fastq a_read1.fastq"))
      system(paste("mv t_read2.fastq a_read2.fastq"))
      system(paste("mv t_read3.fastq a_read3.fastq"))
      
      #####################################################
      #While 2. Run Spades on samples
      #####################################################
      #Run SPADES on sample
      run.contigs<-spades(spades.path = spades.path, "a_read1.fastq", "a_read2.fastq", "a_read3.fastq")
      if (length(run.contigs) > 0) { run.contigs<-run.contigs[width(run.contigs) >= 40] }
      
      if (length(run.contigs) == 0){ break }
      
      #Writes contigs for cap3
      write.loci<-as.list(as.character(run.contigs))
      write.fasta(sequences = write.loci, names = names(write.loci), 
                  "current_seed.fa", nbchar = 1000000, as.string = T)
    
      reference<-"current_seed.fa"
      min.id<-"0.9"
      new.len<-sum(width(run.contigs))

      #####################################################
      #While 4. Checks for endless repeat assemblies
      #####################################################
        
      #If the file gets too large, its due to repeats
      if (new.len >= width(temp.locus)*4){
        
        #Runs CAP3        
        cap3.contigs<-cap3(run.contigs, min.contig.length = 100)

        #Writes contigs for cap3
        save.align<-append(temp.locus, cap3.contigs)
        write.loci<-as.list(as.character(save.align))
        write.fasta(sequences = write.loci, names = names(write.loci), 
                    "temp_seed.fa", nbchar = 1000000, as.string = T)
        
        #Runs MAFFT to align
        system(paste("mafft --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123",
                     " --thread ", threads, " temp_seed.fa > temp_align.fa", sep = ""))
        
        alignment<-scanFa(FaFile("temp_align.fa", sep = ""))   # loads up fasta file
        names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
        
        #Renames sequences to get rid of _R_
        new.align<-strsplit(as.character(alignment), "")
        mat.align<-lapply(new.align, tolower)
        m.align<-as.matrix(as.DNAbin(mat.align))
        
        #Filters out weirdly divergent sequences
        diff<-pairwise.inf.sites(as.character(m.align), loci.names[j])
        bad.seqs<-names(diff)[which(diff >= 0.40)]
        
        #Removes bad contigs
        if (length(bad.seqs) >= 1){ 
          save.contig<-cap3.contigs[!names(cap3.contigs) %in% bad.seqs]
        } else { save.contig<-cap3.contigs }
        
        if (length(save.contig) == 0){
          seeding<-F
          system("rm temp_seed.fa temp_align.fa")
          next
        }
        
        write.loci<-as.list(as.character(save.contig))
        write.fasta(sequences = write.loci, names = names(write.loci), 
                    "current_seed.fa", nbchar = 1000000, as.string = T)
        min.id<-"0.9"
        
        system("rm temp_seed.fa temp_align.fa")
        #makes sure this doesn't go on forever and ever
        repeat.counter<-repeat.counter+1
        if (repeat.counter >= 5){ 
          print(paste("repeat counter hit 5"))
          seeding<-F 
        }#end if
      }#end length > 30,000 if
      
      print(paste("iteration ", counter, " complete!", sep = ""))
      print(paste("new length: ", new.len, ". Old length: ", prev.len, sep = ""))
      if (new.len == prev.len || counter == 20){ 
        seeding<-F 
        print(paste("locus assembly complete after ", counter, " iterations!", sep = ""))
        min.id<-"0.7"
        break
      } #end if 
    }#end while
  
    #####################################################
    #Step 3. Checks contigs to see if there are multiple resulting assemblies
    ##################################################### 
    if (temp.count == 0) { next }
    
    #Checks if there is data  
    if (no.data == T){
      print("No data was found... trying only reads.")
      
      if (file.exists("a_read1.fastq") == T){ 
        system("cat a_read1.fastq  a_read2.fastq  a_read3.fastq >> combined.fasta")
        } else {
        system("cat read1.fq read2.fq singleton.fq >> combined.fasta")
        system("rm read1.fq read2.fq singleton.fq")
      }
      
      #Pulls out the raw reads to just try and use these
      raw.contigs<-sread(readFastq("combined.fa", withIds=TRUE))  # loads up fasta file
      names(raw.contigs)<-paste("seq_", rep(1:length(raw.contigs)), sep = "")
      #Run CAP3
      contigs<-cap3(raw.contigs, min.contig.length = 50)
      system("rm combined.fa")
      
    } else {
      #Save finsihed genome
      contigs<-scanFa(FaFile("current_seed.fa"))   # loads up fasta file
    }
    
    #Skips if there are none
    if (length(contigs) == 0){ next }
    
    #Adds the previous assembly into the contigs if it exists
    if (length(temp.assem) != 0){
      contigs<-append(contigs, temp.assem)
    }
  
    #Trys to merge contigs if there are more than 1
    if (length(contigs) >= 2){
      temp.contigs<-cap3(contigs, min.contig.length = 80)
      contigs<-temp.contigs
      }#end if

    #####################################################
    #Step 4. Checks contigs to see if they actually match to the locus. 
    ##################################################### 
    
    #Align contigs with reference to see which matches
    save.align<-append(temp.locus, contigs)
    final.loci<-as.list(as.character(save.align))
    write.fasta(sequences = final.loci, names = names(final.loci), 
                "saved_contigs.fa", nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste("mafft --localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123",
                 " --thread ", threads, " saved_contigs.fa > loci_align.fa", sep = ""))
    
    alignment<-scanFa(FaFile("loci_align.fa", sep = ""))   # loads up fasta file
    names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
    
    #Renames sequences to get rid of _R_
    new.align<-strsplit(as.character(alignment), "")
    mat.align<-lapply(new.align, tolower)
    m.align<-as.matrix(as.DNAbin(mat.align))
    
    #Filters out weirdly divergent sequences
    diff<-pairwise.inf.sites(as.character(m.align), loci.names[j])
    bad.seqs<-names(diff)[which(diff >= 0.40)]
    
    #Removes bad contigs
    if (length(bad.seqs) >= 1){ 
      save.contig<-contigs[!names(contigs) %in% bad.seqs]
      } else { save.contig<-contigs }
    
   # if (length(save.contig) >= 4) {
    #  break
     # print("many to stitch")
   # }
    
    if (length(save.contig) >= 2){ 
      
      #Align contigs with reference to see which matches
      final.loci<-as.list(as.character(save.contig))
      write.fasta(sequences = final.loci, names = names(final.loci), 
                  "temp_contigs.fa", nbchar = 1000000, as.string = T)
      
      #Matches samples to loci
      system(paste("mpirun pblat -threads=", threads, " reference.fa", 
                   " temp_contigs.fa -tileSize=8 -minIdentity=60", 
                   " -noHead -out=pslx temp_match.pslx", sep = "")) 
      
      #Need to load in transcriptome for each species and take the matching transcripts to the database
      match.data<-fread("temp_match.pslx", sep = "\t", header = F, stringsAsFactors = FALSE)
      setnames(match.data, headers)
      
      #Remove poor, small matches
      match.data<-match.data[match.data$matches > 60,]
      match.data<-match.data[order(match.data$tStart, decreasing = F),]
      con.names<-unique(match.data$qName)
      
      red.match<-c()
      for (k in 1:length(con.names)){
        temp.match<-match.data[match.data$qName == con.names[k],]
        temp.match$tStart[1]<-min(temp.match$tStart)
        temp.match$tEnd[1]<-max(temp.match$tEnd)+(temp.match$qSize[1]-max(temp.match$qEnd))
        red.match<-rbind(red.match, temp.match[1,])
      }
      match.data<-red.match
      merged.seq<-c()
      for (k in 1:nrow(match.data)){
        n.pad<-0
        contig.seq<-save.contig[names(save.contig) == match.data$qName[k]]
        
        #temp.seq<-unlist(strsplit(as.character(final.loci[names(final.loci) == match.data$qName[k]]), ""))
        if (k >= 2) { n.pad<-match.data$tStart[k]-match.data$tEnd[k-1] }
        if (n.pad < 0) { 
          merged.seq<-c()
          break }
        
        if (match.data$strand[k] == "-"){
          seq.string<-reverseComplement(contig.seq)
        } else { seq.string<-contig.seq }
        
        temp.seq<-append(rep("N", n.pad), unlist(strsplit(as.character(seq.string), "")))
        merged.seq<-append(merged.seq, temp.seq)
      }
      
      if (length(merged.seq) != 0){ 
        save.contig<-DNAStringSet(paste(merged.seq, collapse = "")) 
        m.count<-m.count+1
        system("rm temp_match.pslx temp_contigs.fa")
        } else {
        names(save.contig)<-paste(loci.names[j], "_|_",samples[i], "_para_", rep(1:length(save.contig)), sep = "")
        paralog.contigs<-append(paralog.contigs, save.contig)
        print(paste(samples[i], " complete for locus ", loci.names[j], ". These are paralogs.", sep = ""))
        system("rm -r ref current_seed.fa loci_align.fa reference.fa saved_contigs.fa temp_match.pslx temp_contigs.fa")
        next
      }
    }# length(save.contig) > 2 if 
    
    #####################################################
    #Step 5. Saves and cleans up
    ##################################################### 
    if (length(save.contig) == 0){
      print("Attempt failed, no contigs found.")
      next
    }
    
    #Compares new contig length to old contig, saves bigger contigs
    if (length(temp.assem) != 0) {  
      if (width(save.contig) <= width(temp.assem)) { 
        print("Attempt failed, new contig is shorter than old.")
        } else {
        names(save.contig)<-paste(loci.names[j], "_|_",samples[i], sep = "")
        bigger.contigs<-append(bigger.contigs, save.contig)
        print("bigger contig made. ")
        }
      } #end first if
    
    #Checks to see if its new and saves it separately
    if (length(temp.assem) == 0){
      names(save.contig)<-paste(loci.names[j], "_|_",samples[i], sep = "")
      new.contigs<-append(new.contigs, save.contig)
      print("new contig found.")
    }#end if
    
    #Cleanup
    system("rm -r ref current_seed.fa loci_align.fa reference.fa saved_contigs.fa")
    system("rm a_read1.fastq a_read2.fastq a_read3.fastq")
    print(paste(samples[i], " complete for locus ", loci.names[j], sep = ""))
  } #end j loop
    
  #####################################################
  #Step 6. Saves final contigs and cleans up for this species 
  #####################################################  
  write.loci<-as.list(as.character(probe.matched))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(samples[i], "_originals.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  write.loci<-as.list(as.character(paralog.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(samples[i], "_paralogs.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  write.loci<-as.list(as.character(new.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(samples[i], "_new.fa", sep = ""), nbchar = 1000000, as.string = T)  
  
  write.loci<-as.list(as.character(bigger.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(samples[i], "_extended.fa", sep = ""), nbchar = 1000000, as.string = T) 
  
  newer.contigs<-append(new.contigs, bigger.contigs)
  spp.contigs<-probe.matched[!names(probe.matched) %in% names(newer.contigs)]
  final.contigs<-append(spp.contigs, newer.contigs)
  
  write.loci<-as.list(as.character(final.contigs))
  write.fasta(sequences = write.loci, names = names(write.loci), 
              paste(samples[i], "_final.fa", sep = ""), nbchar = 1000000, as.string = T)
  
  system("rm a_read1.fastq a_read2.fastq a_read3.fastq")
   
} #end i loop
  
########################################################################  
# STEP 2
# Output a file of summary stats
########################################################################

#Sets up data summary
header.data<-c("Sample", "finalLoci", "originalLoci", "newLoci", "extendedLoci", "noParalogs", 
               "N50", "N90", "minLen", "maxLen", "meanLen", "medianLen")   
prelim.data<-data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
setnames(prelim.data, header.data)
prelim.data[, Sample:=as.character(samples)]
merge.contigs<-DNAStringSet()

#Cycles through each assembly run and assesses each
for (i in 1:length(samples)){
  
  #Gets raw contigs to count them
  setwd(paste(work.dir, "/", proc.dir, "/", samples[i], "/assembled_loci", sep = ""))
  
  #Gets length of raw contigs
  og<-scanFa(FaFile(paste(samples[i], "_originals.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("originalLoci", header.data), value = length(og) )
  
  #Gets length of raw contigs
  new<-scanFa(FaFile(paste(samples[i], "_new.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("newLoci", header.data), value = length(new) )
  
  #Gets length of raw contigs
  extend<-scanFa(FaFile(paste(samples[i], "_extended.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("extendedLoci", header.data), value = length(extend) )
  
  #Gets length of raw contigs
  para<-scanFa(FaFile(paste(samples[i], "_paralogs.fa", sep = "")))
  temp.para<-unique(gsub("_par.*", "", names(para)))
  set(prelim.data, i = match(samples[i], samples), j = match("noParalogs", header.data), value = length(temp.para) )
  
  #Gets length of raw contigs
  fin<-scanFa(FaFile(paste(samples[i], "_final.fa", sep = "")))
  set(prelim.data, i = match(samples[i], samples), j = match("finalLoci", header.data), value = length(fin) )
  
  #Summarizes basic data
  len.sorted <- rev(sort(width(fin)))
  N50 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.5][1]
  N90 <- len.sorted[cumsum(len.sorted) >= sum(len.sorted)*0.9][1]
  
  set(prelim.data, i =  match(samples[i], samples), j = match("N50", header.data), value = N50 )
  set(prelim.data, i =  match(samples[i], samples), j = match("N90", header.data), value = N90 )
  set(prelim.data, i =  match(samples[i], samples), j = match("minLen", header.data), value = min(width(fin)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("maxLen", header.data), value = max(width(fin)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("meanLen", header.data), value = mean(width(fin)) )
  set(prelim.data, i =  match(samples[i], samples), j = match("medianLen", header.data), value = median(width(fin)) )
  
  merge.contigs<-append(merge.contigs, fin)

}#End loop for things

#Saves combined, final dataset
setwd(work.dir)
write.csv(prelim.data, file = "prealign_extended_sample_assessment.csv")

final.loci<-as.list(as.character(merge.contigs))
write.fasta(sequences = final.loci, names = names(final.loci), "extended_all_species_contigs.fa", nbchar = 1000000, as.string = T)


#### END SCRIPT 
# New files will be located in each species folder within 'Processed_Samples', in the folder assembled_loci

