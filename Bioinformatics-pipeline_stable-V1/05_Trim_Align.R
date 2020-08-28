library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(data.table)

options(stringsAsFactors = FALSE)
#options(warn=2) #for debugging warnings in loops

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

threads<-8
resume<-"yes" #If it crashes it will not redo things already made

###### Parameter setup #########
min.taxa = 4 #min number of taxa to keep an alignment
min.aln = 80 #min size for an alignment
min.len = 60 #min size for an individual sample
edge.trim = 0.51 ### proportion of taxa needed to keep column on ends
min.cov = 0.30 #min proportion of sample that must overlap with consensus

#If you want to remove any taxa
#taxa.remove<-c("Nanorana_parkeri_genome")
taxa.remove = c()

#General directory meanings
work.dir = "/Your/Project/Directory"
align.dir<-"/Your/Project/Directory/Alignments"
probe.file<-"/Your/Project/Directory/Probe_File.fa"
legacy.file<-"/Your/Project/Directory/Final_Legacy_Consensus.fa"
uce.file<-"/Your/Project/Directory/Hutter_uce5k_loci.fa"

# #Example directory for KU cluster
# work.dir = "/home/c111h652/scratch/Mantellidae_All"
# align.dir = "/home/c111h652/scratch/Mantellidae_All/Alignments"
# uce.file = "/home/c111h652/scratch/Mantellidae_All/Hutter_uce5k_loci.fa"
# probe.file = "/home/c111h652/scratch/Mantellidae_All/Master_Ranoidea_All-Markers_Apr21-2019.fa"
# legacy.file = "/home/c111h652/scratch/Mantellidae_All/Final_Legacy_Consensus_Aug4.fa"

###############################################################################
###############################################################################
######################           FUNCTIONS            #########################
###############################################################################
###############################################################################

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

find.orf<-function(input.seq, codons = F, min.size = 80){
  #Sets up data
  # input.seq<-trimmed[j]
  codon.table<-data.frame(Start = rep(0,6), End = rep(0,6), Frame = c("F1", "F2", "F3", "R1", "R2", "R3"))
  for.seq<-as.character(input.seq)
  
  #Gets codon stuff
  TAA<-matchPattern("TAA", for.seq)
  TGA<-matchPattern("TGA", for.seq)
  TAG<-matchPattern("TAG", for.seq)
  
  #Forward Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F1")
    codon.table<-rbind(codon.table, temp.table)
  } 
  
  #Forward Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Forward Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "F3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "F3")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Sets up data
  rev.seq<-as.character(reverseComplement(input.seq))
  
  #Gets codon stuff
  TAA<-matchPattern("TAA", rev.seq)
  TGA<-matchPattern("TGA", rev.seq)
  TAG<-matchPattern("TAG", rev.seq)
  
  #Rev Frame 1
  result1<-TAA[(TAA@ranges@start+2) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+2) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+2) %% 3 == 0]    
  
  starts<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  ends<-c(result1@ranges@start+2, result2@ranges@start+2, result3@ranges@start+2)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R1",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R1")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 2
  result1<-TAA[(TAA@ranges@start+1) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start+1) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start+1) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-1, result2@ranges@start-1, result3@ranges@start-1)
  ends<-c(result1@ranges@start+1, result2@ranges@start+1, result3@ranges@start+1)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R2",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R2")
    codon.table<-rbind(codon.table, temp.table)
  }
  
  #Rev Frame 3
  result1<-TAA[(TAA@ranges@start) %% 3 == 0]   
  result2<-TGA[(TGA@ranges@start) %% 3 == 0]    
  result3<-TAG[(TAG@ranges@start) %% 3 == 0]    
  
  starts<-c(result1@ranges@start-2, result2@ranges@start-2, result3@ranges@start-2)
  ends<-c(result1@ranges@start, result2@ranges@start, result3@ranges@start)
  if (length(starts) != 0){
    codon.table<-codon.table[codon.table$Frame != "R3",]
    temp.table<-data.frame(Start = starts, End = ends, Frame = "R3")
    codon.table<-rbind(codon.table, temp.table)
  } #end if 
  
  if (codons == T) { return(codon.table) }
  
  if (codons == F) {  
    frames<-unique(codon.table$Frame)
    orf.frame<-data.frame()
    for (x in 1:length(frames)){
      temp.codon<-codon.table[codon.table$Frame == frames[x],]
      temp.codon<-temp.codon[order(temp.codon$Start),]
      
      if (temp.codon$Start[1] == 0){
        temp.start<-as.numeric(gsub("F|R", "", temp.codon$Frame))
        add.frame<-data.frame(FrameStart = temp.start, FrameEnd = width(input.seq), 
                              Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
        orf.frame<-rbind(orf.frame, add.frame)
        next
      }
      #Goes through each of the given directions codons and converts to frame ranges
      temp.frame<-data.frame()
      for (y in 1:(nrow(temp.codon)+1)){
        #First y the start is 1, otherwise take from previous end
        if (y == 1){ frame.start<-as.numeric(gsub("F|R", "", temp.codon$Frame[y])) } else { frame.start<-temp.frame$FrameEnd[y-1]+4 }
        
        #Gets end by subtracting from the codon start
        frame.end<-temp.codon$Start[y]-1
        temp.frame<-rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
      } # end y loop
      
      temp.frame$FrameEnd[nrow(temp.frame)]<-width(input.seq)
      
      #Adds all the data together
      add.frame<-cbind(temp.frame, Size = (temp.frame$FrameEnd-temp.frame$FrameStart)+1, Frame = frames[x])
      orf.frame<-rbind(orf.frame, add.frame)
      
    } #end x loop
    
    orf.frame<-orf.frame[orf.frame$Size >= min.size,]
    return(orf.frame)
  } # end else
  
}# END FUNCTION

trim.ends<-function (x, min.n.seq = 4, codon.trim = T){
  #Converts DNAStringSet to something usable
  #x<-trimmed
  #min.n.seq<-4
  new.align<-strsplit(as.character(x), "")
  mat.align<-lapply(new.align, tolower)
  x<-as.matrix(as.DNAbin(mat.align))

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
  #Does stuff
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b) {
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if (max(m) < min.n.seq) { stop("alignment contains less sequences then required") }
  m <- range(which(m >= min.n.seq))
  
  #Forward Frame 2
  if (codon.trim == T){
    if ((m[1]-1) %% 3 == 0){ m[1]<-m[1] }
    if ((m[1]-1) %% 3 == 1){ m[1]<-m[1]+2 }
    if ((m[1]-1) %% 3 == 2){ m[1]<-m[1]+1 }
  }
  
  m <- seq(from = m[1], to = m[2])
  x2 <- as.matrix(x[, m])
  #Converts back
  save.names<-rownames(x2)
  
  #Removes N end gaps
  x3<-as.list(data.frame(t(as.character(x2))))
  for (y in 1:length(x3)){
   #Starts from the beginning and end to fill in end gaps
    for (q in 1:length(x3[[y]])){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
    for (q in length(x3[[y]]):1){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
  }#end x loop
  #Saves final stuff
  temp.align<-lapply(x3, FUN = function(x) paste(x, collapse = ""))
  align.out<-DNAStringSet(unlist(temp.align))
  names(align.out)<-save.names
  return(align.out)
}

run.mafft<-function(unaligned.contigs, add.contigs = NULL, algorithm = "localpair", rev.dir = T, save.name = NULL, threads = 6, delete.files = T){
  
  #save.name<-locus.save.name
  #algorithm = "localpair"
  # unaligned.contigs<-intron.align
  save.contigs<-as.list(as.character(unaligned.contigs))
  if (is.null(save.name) == T) { save.name<-paste(sample(LETTERS, 5, replace = T), collapse = "")}
  if (rev.dir == T){ adjust.direction<-"--adjustdirection" } else { adjust.direction<-"" }
  
  #Adds a sequence into the alignment. Saves much computation.
  if (algorithm == "add"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs), 
                paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)
    
    #Saves to folder to run with mafft
    add.save<-as.list(as.character(add.contigs))
    write.fasta(sequences = add.save, names = names(add.save), 
                "add_sequences.fa", nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste("mafft --",algorithm, " add_sequences.fa ", adjust.direction, " --maxiterate 1000 ", 
                 save.name, ".fa > ", save.name, "_align.fa", sep = ""), ignore.stderr = T)
    
    alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
    unlink("add_sequences.fa")
    
  }#end -add
  
  #Does Regular MAFFT Local Pair
  if (algorithm == "localpair"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs), 
                paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)
    
    #Runs MAFFT to align
    system(paste("mafft --",algorithm, " --maxiterate 1000 ", adjust.direction, " --quiet --op 3 --ep 0.123",
                 " --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa", sep = ""))
    
    alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
  }#end local pair
  
  if (delete.files == T){
    unlink(paste(save.name, "_align.fa", sep = ""))
    return(alignment)
  } else { return(alignment) }
}#function end

pairwise.inf.sites<-function(x, y) {
  #Alignment should be DNAStringSet
  # x<-m.align
  # y<-"Reference_Locus"
  temp.align<-strsplit(as.character(x), "")
  mat.align<-lapply(temp.align, tolower)
  m.align<-as.matrix(as.DNAbin(mat.align))
  
  #Filters out weirdly divergent sequences
  new.align<-as.character(m.align)
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

slice.trim<-function(input.align, slice.size.bp = 100, threshold = 0.45){
  
  #makes consensus sequence for comparison
  #input.align<-trimal.align
  input.con<-make.consensus(input.align, method = "majority")
  names(input.con)<-"Reference_Locus"
  
  comb.align<-append(input.align, input.con)
  
  #Gets slice information ready
  slice.no<-ceiling(max(width(input.align))/slice.size.bp)
  slice.start<-1
  slice.end<-slice.size.bp
  
  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){ 
    slice.end<-max(width(input.align))
  }#end if check
  output.align<-DNAStringSet()
  for (x in 1:slice.no){
    
    #Slice alignment into number of slices 
    sliced.align<-subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences 
    bad.align<-pairwise.inf.sites(sliced.align, "Reference_Locus")
    #Remove bad sequence chunks
    rem.seqs<-bad.align[bad.align >= threshold]
    good.align<-sliced.align[!names(sliced.align) %in% names(rem.seqs)]
    #Makes replacement gap seqs for the bad ones
    blank.align<-DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if
    
    #Saves the slices and cats
    save.slice<-append(good.align, blank.align)
    save.slice<-save.slice[order(names(save.slice))]
    save.names<-names(save.slice)
    output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align)<-save.names
    
    #Gets new start and stop
    slice.start<-slice.start+100
    slice.end<-slice.end+100
    #checks to see if its out of bounds
    if (slice.end > max(width(input.align))){ 
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){ break } else {
        save.slice<-subseq(comb.align, start = slice.start, end = slice.end)
        save.slice<-save.slice[order(names(save.slice))]
        save.names<-names(save.slice)
        output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
        names(output.align)<-save.names
        break
        }
    }#end if
  }#end x loop
  
  #Removes reference
  output.align<-output.align[names(output.align) != "Reference_Locus"]  
  #removes gap only taxa
  str.splitted<-strsplit(as.character(output.align), "")
  x.align<-as.matrix(as.DNAbin(str.splitted) )
  len.temp<-as.character(as.list(x.align))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= 20]  
  return.align<-output.align[!names(output.align) %in% names(spp.rem)]
  return(return.align)      
}#end FUNCTION

run.trimal<-function(input.align, method = "auto"){
  #input.align<-rem.align
  #Finds probes that match to two or more contigs
  save.rownames<-names(input.align)
  write.align<-as.list(as.character(input.align))
  write.fasta(sequences = write.align, names = names(write.align), 
              paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".fa"), nbchar = 1000000, as.string = T)
  
  input.file<-paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".fa")
  
  system(paste0("trimal -in ", input.file, " -out ", input.file, "-tm -automated1"))
  system(paste0("rm ", input.file))
  
  if (file.exists(paste0(input.file, "-tm")) == F) { 
    system(paste0("rm ", input.file))
    print(paste0("deleted. Not enough overlapping data in alignment.") )
    return(DNAStringSet())
  } else { system(paste0("mv ", input.file, "-tm ", input.file)) }
  
  out.align<-scanFa(FaFile(input.file)) 
  
  #Fixes any terrible NA names introduced by trimal
  new.names<-c()
  for (j in 1:length(names(out.align))){ 
    new.names[j]<-save.rownames[grep(pattern = names(out.align)[j], x = save.rownames)]
  }
  
  temp<-names(out.align)[is.na(names(out.align)) == T]
  if (length(temp) > 0){ stop("there are NAs in the names") }
  names(out.align)<-new.names
  unlink(input.file)
  return(out.align)
}#end function

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


################################################################################################
################################################################################################
############## Step 1. Trim the exon-intron alignments combined     ############################
################################################################################################
################################################################################################
################################################################################################

trim.dir<-"all-markers_trimmed"
dir.create(paste0(align.dir, "/", trim.dir))
setwd(paste0(align.dir, "/all-markers_untrimmed"))
locus.names<-list.files(".")

if (resume == "yes"){
  done<-list.files(paste0(align.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
 
  ##############
  #STEP 1: Setup steps
  ##############
  #Reads in files
  setwd(paste0(align.dir, "/all-markers_untrimmed"))
  align<-DNAStringSet(readAAMultipleAlignment(file = paste0(align.dir, "/all-markers_untrimmed/", locus.names[i]), format = "phylip"))

  #Remove gap only alignments
  gap.align<-strsplit(as.character(align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  rem.align<-align[!names(align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(rem.align) <= as.numeric(min.taxa)){ next  }
  #removes too short loci
  if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 2: Runs trimming programs
  ##############
  #A. Runs Trimal
  setwd(paste0(align.dir, "/", trim.dir))
  trimal.align<-run.trimal(rem.align)

  #removes too short loci
  if (max(width(trimal.align)) <= as.numeric(min.aln)){ next }
  
  #B Slice up alignment
  if (max(width(trimal.align)) >= 100){
    red.align<-slice.trim(trimal.align, slice.size.bp = 80, threshold = 0.40)
    if (length(red.align) == 0){ next }
  }else {red.align = trimal.align }
  
  #C. Trims the alignment ends
  edge.align<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){ next }
  
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(edge.align, method = "majority")
  #Removes the edge gaps
  ok.seq<-gsub("\\+|-", "", as.character(con.seq))
  #removes loci with too few taxa
  if ( nchar(ok.seq) <= as.numeric(min.aln)){ next  }

  ##############
  #STEP 3: Cleanup and save
  ##############
  #removes loci with too few taxa
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste0(locus.names[i], " deleted. Too few taxa after trimming.") )
    next
  }
  
  write.temp<-strsplit(as.character(edge.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)

}
  
################################################################################################
################################################################################################
############## Step 2. Verify and retrimm exons. Output exon table  ############################
################################################################################################
################################################################################################
################################################################################################

#Sets up new directory for this stuff
trim.dir<-"exon-only_trimmed"
dir.create(paste0(align.dir, "/", trim.dir))
prot.dir<-"exon-protein_trimmed"
dir.create(paste0(align.dir, "/", prot.dir))
setwd(paste0(align.dir, "/exon-only_untrimmed"))
locus.names<-list.files(".")

#Checks if to resume or not
if (resume == "yes"){
  done<-list.files(paste0(align.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Load in data and basic filtering
  ##############
  #Skip if a UCE
  if (length(grep("uce", locus.names[i])) == 1){ next }
  
  #Reads in files
  setwd(paste0(align.dir, "/exon-only_untrimmed"))
  align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))

  #Remove gap only sequences
  gap.align<-strsplit(as.character(align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  red.align<-align[!names(align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(red.align) <= as.numeric(min.taxa)){ next }
  #removes too short loci
  if (max(width(red.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 2: Trims the alignment 
  ##############  
  trimmed<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = T)
  save.names<-names(trimmed)
  
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(trimmed, method = "majority")
  
  #Removes the edge gaps
  ref.aligned<-as.character(con.seq)
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]

  #Finds weird gaps to fix
  temp.gaps<-as.numeric(1)
  for (k in 1:length(not.gaps)-1){ temp.gaps<-append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
  temp.gaps<-temp.gaps-1
  names(temp.gaps)<-not.gaps
  gap.spots<-temp.gaps[temp.gaps %% 3 != 0]
  
  del.col<-c()
  if (length(gap.spots) != 0){
    #Loop through each potential bad thing and fix
    for (k in 1:length(gap.spots)){
      del.col<-append(del.col, (as.numeric(names(gap.spots[k]))-gap.spots[k]):(as.numeric(names(gap.spots[k]))-1))
    }
  }#end if
  
  ##############
  #STEP 3: Fixes large gaps at ends of alignment
  ##############  
  #Looks for gaps that clearly wrong and not 3 BP 
  new.align<-strsplit(as.character(trimmed), "")
  x<-as.matrix(as.DNAbin(new.align))
  
  rem.n<-c()
  for (k in 1:ncol(x)){
    gaps<-table(as.character(x[,k]))
    per.gaps<-gaps[names(gaps) == "-"]/nrow(x)
    
    if (length(per.gaps) == 0){ next }
    
    #Records column when the gaps exceed this percentage
    if (per.gaps >= 0.75){ del.col<-append(del.col, k) }
  
    #Removes gap columns only consisting of Ns 
    n.gaps<-gaps[names(gaps) != "-"]
    if (length(n.gaps) == 1){
      if (names(n.gaps) == "n"){ rem.n<-append(rem.n, k)}
    }
    
  }#end k loop

  #combines columns to be deleted
  fin.del<-c(rem.n, del.col)
  if (length(fin.del) != 0){ x<-x[,-fin.del] }
  #Removes bad columsn and coverts alignment back to DNASTringSet
  char.align<-as.list(data.frame(t(as.character(x))))
  temp.align<-lapply(char.align, FUN = function(x) paste(x, collapse = ""))
  trimmed<-DNAStringSet(unlist(temp.align))
  names(trimmed)<-save.names

  ##############
  #STEP 4: Gathers table of best and longest stop codon free frames for each seq
  ##############  
  #Checks to make sure the codon position is correct
  save.frame<-data.frame()
  save.all<-data.frame()
  for (j in 1:length(trimmed)){
    #Finds open reading frames
    temp.codon<-find.orf(trimmed[j], codons = F, min.size = 80 )
    if(nrow(temp.codon) == 0){
      samp.frame<-cbind(Sample = names(trimmed[j]), FrameStart = 0, FrameEnd = 0, Size = 0, sppSize = 0, Frame = "0")
      save.frame<-rbind(save.frame, samp.frame)
      next
    }
    
    all.frame<-temp.codon[temp.codon$Size >= max(temp.codon$Size) * .70,]
    big.frame<-temp.codon[temp.codon$Size == max(temp.codon$Size),]

    if (nrow(big.frame) >= 2){
      #Picks the best from this order of things
      temp.stop<-big.frame[big.frame$Frame == "F1",]
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R1",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F2",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R2",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F3",] }
      if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R3",] }
      big.frame<-temp.stop
    }

    # #Saves teh data
    samp.frame<-cbind(Sample = names(trimmed[j]), big.frame)
    temp.size<-unlist(strsplit(as.character(trimmed[j]), ""), use.names = F)
    
    #Starts from the beginning and end to fill in end gaps
    sub.size<-0
    for (q in 1:length(temp.size)){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
    for (q in length(temp.size):1){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
    
    #Saves final data
    samp.frame<-cbind(samp.frame, sppSize = length(temp.size)-sub.size)
    save.frame<-rbind(save.frame, samp.frame)
    all.frame<-cbind(Sample = names(trimmed[j]), all.frame)
    all.frame<-cbind(all.frame, sppSize = length(temp.size)-sub.size)
    save.all<-rbind(save.all, all.frame)
  }#end j loop
  
  #Moves on if there are no frames found. Saves to Anon folder? 
  if (unique(save.frame$Frame)[1] == "0"){
    print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
    next
  }
  
  ##############
  #STEP 5: Uses previous data to find a consistent frame
  ##############  
  #Looks at the overall data rather than the best indiv data to find a consistent frame
  temp.all<-save.all
  frame.names<-unique(temp.all$Frame)
  #Goes through the equally good frames and reduces to frames with the same range
  very.best<-data.frame()
  for (k in 1:length(frame.names)){
    temp.best<-temp.all[temp.all$Frame == frame.names[k],]
    starts<-table(temp.best$FrameStart)[table(temp.best$FrameStart) == max(table(temp.best$FrameStart))]
    ends<-table(temp.best$FrameEnd)[table(temp.best$FrameEnd) == max(table(temp.best$FrameEnd))]
    
    #Removes duplicates
    starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
    ends<-ends[as.numeric(names(ends)) == min(as.numeric(names(ends)))]
    
    super.best<-temp.best[temp.best$FrameStart == as.numeric(names(starts)),]
    super.best<-super.best[super.best$FrameEnd == as.numeric(names(ends)),]
    very.best<-rbind(very.best, super.best)
  }#end k loop
  
  #Moves on if there are no frames found. Saves to Anon folder? 
  if (nrow(very.best) == 0){
    print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
    next
  }

  ##############
  #STEP 6: Selects the best frame
  ##############  
  #Picks out the best frame
  best.frame<-table(very.best$Frame)[table(very.best$Frame) == max(table(very.best$Frame))]
  
  #If there are multiple good frames pick the biggest
  if (length(best.frame) != 1){ 
    temp.fix<-very.best[very.best$Frame %in% names(best.frame),]
    bigger<-temp.fix[temp.fix$Size == max(temp.fix$Size),]
    best.frame<-table(bigger$Frame)[table(bigger$Frame) == max(table(bigger$Frame))]
  }#end if
  
  #If they are same size just pick from this order
  if (length(best.frame) != 1){ 
    #Picks the best from this order of things
    temp.stop<-best.frame[names(best.frame) == "F1"]
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R1"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F2"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R2"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F3"] }
    if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R3"] }
    best.frame<-temp.stop
    }
    
  #Checks the remaining ORF size
  temp.size<-save.all[save.all$Frame == names(best.frame),]
  samp.spp<-temp.size[temp.size$sppSize == max(temp.size$sppSize),]
  samp.seq<-trimmed[names(trimmed) == samp.spp$Sample[1]]
  samp.size<-nchar(gsub("-", "", as.character(samp.seq)))
  
  if (mean(temp.size$Size) <= samp.size*.5){ print(paste(locus.names[i], " was small.", sep = "")) }
  
  if (mean(temp.size$Size) <= samp.size*.25){ 
    print(paste(locus.names[i], " sucked. Too few sequence left.", sep = ""))
    next
  }
  
  if (best.frame <= length(trimmed) * .5){ 
    print(paste(locus.names[i], " sucked. No cosistent frame.", sep = ""))
    next
  }
  
  #Reverses if it needs to
  if (length(grep("R", names(best.frame))) != 0){
    new.align<-reverseComplement(trimmed)
  }else { new.align<-trimmed }
  
  ##############
  #STEP 7: Gets start and stop coordinates for each sequence and find best alignment
  ##############  
  #Gets trimming locations
  frame.ranges<-save.all[save.all$Frame == names(best.frame),]
  
  #Gets potential starts and ends
  starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
  ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
  starts<-starts[starts == max(starts)]
  ends<-ends[ends == max(ends)]
  
  if (length(starts) != 1 || length(ends) != 1){ 
    frame.ranges<-save.all[save.all$Frame == names(best.frame),]
    starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
    ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
    starts<-starts[starts == max(starts)]
    ends<-ends[ends == max(ends)]
  }
  
  if (length(starts) != 1 || length(ends) != 1){ 
    starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
    ends<-ends[as.numeric(names(ends)) == max(as.numeric(names(ends)))]
  }
  
  ###################
  #STEP 8: Makes sure entire alignment is a multiple of 3
  ###################
  anu.start<-as.numeric(names(starts))
  new.end<-as.numeric(names(ends))
  new.len<-new.end-(anu.start-1)
  
  #Gets a new end to keep in multiples of 3 for proteins
  if (length(new.len[which(new.len %%3==0)]) == 0) {
    anu.end<-new.end-1
  } else { anu.end<-new.end }
  
  new.len<-anu.end-(anu.start-1)
  if (length(new.len[which(new.len %%3==0)]) == 0) {
    anu.end<-new.end-2
  } else { anu.end<-anu.end }
  
  #Trims sequence with new coords
  done.seq<-subseq(start = anu.start, end = anu.end, x = new.align)

  ###################
  #STEP 9: Trim out odd start/end bases
  ###################
  codon.seq<-DNAStringSet()
  for (k in 1:length(done.seq)){
    ref.aligned<-as.character(done.seq[k])
    
    #Chcecks at beginning of sequence
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    
    if (length(not.gaps) <= min.len){ next }
    
    #Checks if its odd, delete 1 base
    if ( (not.gaps[1]-1) %%3 == 2){ substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-" }
    
    #Deletes 2 bases its off by
    if ( (not.gaps[1]-1) %%3 == 1){
      substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-"
      substr(ref.aligned, not.gaps[1]+1, not.gaps[1]+1)<-"-"
    }#end if
    
    #checks for end of sequence
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    #counts characters
    char.len<-(not.gaps[length(not.gaps)]-not.gaps[1])+1
    end.pos<-not.gaps[length(not.gaps)]
    #removes odd characters at ends
    if ( char.len %%3 == 1){ substr(ref.aligned, end.pos, end.pos)<-"-" }
    
    if ( char.len %%3 == 2){ 
      substr(ref.aligned, end.pos-1, end.pos-1)<-"-"
      substr(ref.aligned, end.pos, end.pos)<-"-"
    } #end if 
    
    #Saves final seqs
    save.seq<-DNAStringSet(ref.aligned)
    names(save.seq)<-names(done.seq[k])
    codon.seq<-append(codon.seq, save.seq)
    
  } #END K  
  
  ###################
  #STEP 10: Change stop codons to N
  ###################
  #Finds stop codons to replace
  n.seq<-DNAStringSet(gsub("-", "N", as.character(codon.seq)))
  stop.seq<-DNAStringSet()
  for (k in 1:length(n.seq)){
    stop.data<-find.orf(n.seq[k], codon = T, min.size = 80)
    stop.data<-stop.data[stop.data$Frame == "F1",]
    
    #Skips if there are more than 3 stop codons
    if (nrow(stop.data) >= 3){ next }
    
    if (stop.data$Start[1] == 0){ stop.seq<-append(stop.seq, n.seq[k]) } else {
      #Goes through each codon
      ref.aligned<-as.character(n.seq[k])
      for (y in 1:nrow(stop.data)){      
        #Saves final seqs
        substr(ref.aligned, stop.data$Start[y], stop.data$Start[y])<-"N"
        substr(ref.aligned, stop.data$Start[y]+1, stop.data$Start[y]+1)<-"N"
        substr(ref.aligned, stop.data$Start[y]+2, stop.data$Start[y]+2)<-"N"
      }#end Y LOOP
      
      #Saves final data
      save.seq<-DNAStringSet(ref.aligned)
      names(save.seq)<-names(n.seq[k])
      stop.seq<-append(stop.seq, save.seq)
    }#end if state
   }# END K loop
  
  if (length(stop.seq) <= min.taxa){ next }
  
  ###################
  #FINAL STEP: Save everything after some final spp and length filtering
  ###################
  #Removes sequences that are less than a certain coverage
  t.align<-strsplit(as.character(stop.seq), "")
  len.loci<-lapply(t.align, function (x) x[x != "N"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  
  if (length(spp.rem) > 0){ 
    red.align<-t.align[!names(t.align) %in% unique(names(spp.rem))]
  } else { red.align<-t.align }
  
  #Removes if all removed
  if (length(red.align) == 0){ 
    print(paste0("not enough alignment remains for ", locus.names[i]))
    next
  }
  
  #writes alignment
  mat.align<-lapply(red.align, tolower)
  write.align<-as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(write.align, file=paste0(align.dir, "/", trim.dir, "/", locus.names[i]), interleave = F)
  
  #Saves protein sequence version
#  trans.prot<-translate(stop.seq, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
#  mat.align<-lapply(trans.prot, tolower)
#  write.align<-as.matrix(as.AAbin(trans.prot))
#  write.phy(write.align, file=paste0(align.dir, "/", prot.dir, "/", locus.names[i]), interleave = F)
  
  #CHECKS FOR STOP CODNS
  #stop.align<-strsplit(as.character(trans.prot), "")
  #stop.loci<-lapply(stop.align, function (x) x[x == "*"])
  #stop.len<-unlist(lapply(stop.loci, function (x) length(x)))
  #stop.pres<-stop.len[stop.len >= 1]  
  
  #if (length(stop.pres) != 0){ stop("STOP CODON") }
  #IQTREE TEST
 # system(paste0("iqtree -s ", align.dir, "/", trim.dir, "/", locus.names[i], 
  #             " -nt ", threads, " -m MFP -st CODON -rcluster 10 -msub nuclear"))
  
  
}#end i loop
  

################################################################################################
################################################################################################
############## Step 3. Save UCE and bad exon loci separately        ############################
################################################################################################
################################################################################################
################################################################################################

#Blast two probe set files together to get uce loci to keep
setwd(align.dir)

#Make blast database for the probe loci
system(paste("makeblastdb -in ", probe.file, " -parse_seqids -dbtype nucl ",
             " -out probe_blast_db", sep = ""))

#Matches samples to loci
system(paste("blastn -task dc-megablast -db probe_blast_db",
             " -query ", uce.file, " -out uce_match.txt", 
             " -outfmt 6 -num_threads ", threads, sep = ""))

#headers for the blast db
headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
           "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

#Load in matches
match.data<-fread("uce_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
setnames(match.data, headers)
match.data<-match.data[match.data$evalue <= 0.05,]
uce.names<-unique(match.data$tName)
uce.names<-paste(uce.names, ".phy", sep = "")
system("rm probe_blast* uce_match.txt")

#Sets up new directory for this stuff
trim.dir<-"uce_trimmed"
dir.create(paste(align.dir, "/", trim.dir, sep = ""))
setwd(paste(align.dir, "/all-markers_untrimmed", sep = ""))
file.names<-list.files(".")
locus.names<-uce.names[uce.names %in% file.names]

#Checks if to resume or not
if (resume == "yes"){
  done<-list.files(paste0(align.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Basic steps
  ##############
  #Reads in files
  setwd(paste(align.dir, "/all-markers_untrimmed", sep = ""))
  align<-readAAMultipleAlignment(file = paste(align.dir, "/all-markers_untrimmed/", locus.names[i], sep =""), format = "phylip")
  align<-DNAStringSet(align)
  
  #Remove gap only alignments
  gap.align<-strsplit(as.character(align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  rem.align<-align[!names(align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(rem.align) <= as.numeric(min.taxa)){ next }
  
  #removes too short loci
  if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 2: Runs trimming programs
  ##############
  #A. Runs Trimal
  trimal.align<-run.trimal(rem.align)
  
  if (max(width(trimal.align)) <= as.numeric(min.aln)){ next }
  
  #B Slice up alignment
  red.align<-slice.trim(trimal.align, slice.size.bp = 100, threshold = 0.40)
  
  #C. Trims the alignment ends
  edge.align<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(edge.align, method = "majority")
  #Removes the edge gaps
  ok.seq<-gsub("\\+|-", "", as.character(con.seq))
  
  #removes loci with too few taxa
  if ( nchar(ok.seq) <= as.numeric(min.aln)){ next  }
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 3: Cleanup and save
  ##############
  #removes loci with too few taxa
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = "") )
    next
  }
  #string splitting
  write.temp<-strsplit(as.character(edge.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(align.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}


################################################################################################
################################################################################################
############## Step 4. Save the intron regions along separately     ############################
################################################################################################
################################################################################################
################################################################################################

### GET EXON_INTRON FILES NOT USED FOR EXON 

#Sets up new directory for this stuff
setwd(align.dir)
trim.dir<-"intron-only_trimmed"
dir.create(paste(align.dir, "/", trim.dir, sep = ""))
notrim.dir<-"intron-only_untrimmed"
dir.create(paste(align.dir, "/", notrim.dir, sep = ""))
setwd(paste(align.dir, "/exon-only_trimmed", sep = ""))
exon.names<-list.files(".")

#Checks nad removes some missing from others
setwd(paste(align.dir, "/all-markers_untrimmed", sep = ""))
all.names<-list.files(".")
locus.names<-all.names[all.names %in% exon.names]

#Checks if to resume or not
if (resume == "yes"){
  done<-list.files(paste0(align.dir, "/", trim.dir))
  locus.names<-locus.names[!locus.names %in% done]
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Set up and match loci
  ##############
  #Reads in files
  setwd(paste(align.dir, "/exon-only_trimmed", sep = ""))
  align<-readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  align<-DNAStringSet(align)
  
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(align, method = "majority")
  #Removes the edge gaps
  ok.seq<-DNAStringSet(gsub("\\+|-", "", as.character(con.seq)))
  names(ok.seq)<-paste("Reference_Locus")
  
  setwd(paste(align.dir, "/all-markers_untrimmed", sep = ""))
  intron.align<-DNAStringSet(readAAMultipleAlignment(file = locus.names[i], format = "phylip"))

  ##############
  #STEP 2: Runs MAFFT to add
  ##############
  setwd(paste(align.dir, "/", trim.dir, sep = ""))
  alignment<-run.mafft(unaligned.contigs = intron.align, add.contigs = ok.seq, rev.dir = T,
                       algorithm = "add", save.name = gsub(".phy", "",locus.names[i]), delete.files = T)
  
  #Aligns and then reverses back to correction orientation
  reversed<-names(alignment)[grep(pattern = "_R_", names(alignment))]
  if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment<-reverseComplement(alignment) }
  names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  
  #Gets the divergence to make sure not crazy
  diff<-pairwise.inf.sites(alignment, "Reference_Locus")
  bad.seqs<-names(diff)[which(diff >= 0.5)]
  rem.align<-alignment[!names(alignment) %in% bad.seqs]
  
  # Moves onto next loop in there are no good sequences
  if (length(rem.align) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " had too few taxa", sep = ""))
    next }
  
  #removes too short loci
  #if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 3: Removes exon from the intron part
  ##############
  #Removes the edge gaps
  ref.aligned<-as.character(alignment['Reference_Locus'])
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start<-min(not.gaps)
  ref.finish<-max(not.gaps)
  
  #Finds weird gaps to fix
  temp.gaps<-as.numeric(1)
  for (k in 1:length(not.gaps)-1){ temp.gaps<-append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
  temp.gaps<-temp.gaps-1
  names(temp.gaps)<-not.gaps
  bad.gaps<-which(temp.gaps >= 30)
  front.gaps<-bad.gaps[bad.gaps <= length(not.gaps) *0.10]
  end.gaps<-bad.gaps[bad.gaps >= length(not.gaps) *0.90]

  #Fix big gaps if there are any
  if (length(front.gaps) != 0){ 
    temp.change<-(max(as.numeric(names(front.gaps))-ref.start))-(max(front.gaps)-1)
    ref.start<-ref.start+temp.change
  }#end gap if
  
  #Fix big gaps if there are any
  if (length(end.gaps) != 0){ 
    add.bp<-length(temp.gaps)-min(end.gaps)
    #add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
    min.gaps<-temp.gaps[min(end.gaps)]
    temp.change<-as.numeric(names(min.gaps))-as.numeric(min.gaps)
    ref.finish<-temp.change+add.bp
  }#end gap if
  
  #Cuts out the intron pieces
  intron.left<-subseq(alignment, 1, ref.start-1)
  intron.right<-subseq(alignment, ref.finish+1, width(alignment))
  save.names<-names(alignment)
  
  #Merges the alignments
  intron.align<-DNAStringSet(paste0(as.character(intron.left), as.character(intron.right)))
  names(intron.align)<-save.names
  intron.align<-intron.align[names(intron.align) != "Reference_Locus"]
  
  #Remove gap only alignments
  gap.align<-strsplit(as.character(intron.align), "")
  gap.count<-unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem<-gap.count[gap.count <= as.numeric(min.aln)]
  
  #Remove the taxa that need to go
  rem.taxa<-c(names(gap.rem), taxa.remove)
  rem.align<-intron.align[!names(intron.align) %in% rem.taxa]
  
  #removes loci with too few taxa
  if ( length(rem.align) <= as.numeric(min.taxa)){ next }
  
  #removes too short loci
  if (max(width(rem.align)) <= as.numeric(min.aln)){ next }
  
  #Saves a copy of not modified intron only
  intron.notrim<-strsplit(as.character(rem.align), "")
  save.intron<-as.matrix(as.DNAbin(intron.notrim) )
  write.phy(save.intron, file= paste0(align.dir, "/", notrim.dir, "/", locus.names[i]), interleave = F)
  
  ##############
  #STEP 4: Runs trimming programs
  ##############
  #A. Runs Trimal
  trimal.align<-run.trimal(rem.align)
  
  if (length(trimal.align) == 0){ next }
  
  #B Slice up alignment
  red.align<-slice.trim(trimal.align, slice.size.bp = 100, threshold = 0.5)
  
  if (length(red.align) == 0){ next }
  
  #C. Trims the alignment ends
  edge.align<-trim.ends(red.align, min.n.seq = ceiling(length(red.align) * edge.trim), codon.trim = F)
  #Gets consensus seq for trimming more
  con.seq<-make.consensus(edge.align, method = "majority")
  #Removes the edge gaps
  ok.seq<-gsub("\\+|-", "", as.character(con.seq))
  
  #removes loci with too few taxa
  if ( nchar(ok.seq) <= as.numeric(min.aln)){ next  }
  #removes too short loci
  if (max(width(edge.align)) <= as.numeric(min.aln)){ next }
  
  ##############
  #STEP 5: Cleanup and save
  ##############
  #removes loci with too few taxa
  if (length(names(edge.align)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = "") )
    next
  }
  #string splitting
  write.temp<-strsplit(as.character(edge.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #removes too short loci
  if (ncol(aligned.set) <= as.numeric(min.aln)){ 
    print(paste(locus.names[i], " deleted. Trimmed alignment length below threshold.", sep = "") )
    next 
  }
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(align.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}


################################################################################################
################################################################################################
############## Step 5. Group linked loci and exons from same gene   ############################
################################################################################################
################################################################################################
################################################################################################

#Sets up new directory for this stuff
trim.dir<-"locus-combined"
exon.dir<-"exon-only_trimmed"
dir.create(paste(align.dir, "/", trim.dir, sep = ""))
setwd(paste(align.dir, "/", exon.dir, sep = ""))

#Finds loci from the file names
file.names<-list.files(".")
exon.names<-file.names[grep("-ex", file.names)]
temp.names<-gsub(".*_", "", exon.names)
temp.names<-gsub("-.*", "", temp.names)
temp.names<-gsub(".phy", "", temp.names)
locus.names<-unique(temp.names[duplicated(temp.names)])

#Fixes broken loci
if (length(locus.names) <= 10){
  exon.names<-file.names[grep("-ex", file.names)]
  temp.names<-str_split(exon.names, "-")
  locus.names <- unique(unlist(lapply(temp.names, "[[", 3)))
}

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){
  
  ##############
  #STEP 1: Get the taxa in the alignments for these loci
  ##############
  #Current probe sets with hyphens
  if (length(grep(".phy$", locus.names[i])) >= 1){ next }
  setwd(paste(align.dir, "/", exon.dir, sep = ""))
  locus.files<-exon.names[grep(paste0("-", locus.names[i], "-"), exon.names)]
  locus.files<-locus.files[order(gsub(".*-ex", "", locus.files))]
  
  #Older probe sets with underscore
  if (length(locus.files) == 0){
    locus.files<-exon.names[grep(paste0("_", locus.names[i], "-"), exon.names)]
    locus.files<-locus.files[order(gsub(".*-ex", "", locus.files))]
  } #end if
  
  if (length(grep("-ex2.phy", locus.files)) > 1){ stop("grep FAIL") }
  
  taxa.names<-c()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste(align.dir, "/", exon.dir, "/", locus.files[j], sep =""), format = "phylip")
    taxa.names<-append(taxa.names, rownames(align))
  }
  
  taxa.names<-unique(taxa.names)
  
  ##############
  #STEP 2: Get different loci
  ##############
  
  combined.align<-DNAStringSet()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste(align.dir, "/", exon.dir, "/", locus.files[j], sep =""), format = "phylip")
    align<-DNAStringSet(align)

    add.taxa<-taxa.names[!taxa.names %in% names(align)]
    
    blank.align<-DNAStringSet()
    if (length(add.taxa) != 0){
      for (y in 1:length(add.taxa)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", max(width(align))), collapse = "")) )
      }
      names(blank.align)<-add.taxa
    }#end rem seqs if
    
    #Saves the slices and cats
    new.align<-append(align, blank.align)
    new.align<-new.align[order(names(new.align))]
    save.names<-names(new.align)
    combined.align<-DNAStringSet(paste0(as.character(combined.align), as.character(new.align)))
    names(combined.align)<-save.names
  }#end j loop
  
  ##############
  #STEP 3: Cleanup and save
  ##############

  #string splitting
  write.temp<-strsplit(as.character(combined.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )
  
  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]  
  spp.rem<-append(spp.rem, spp.len[spp.len <= min.len]  )
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }
  
  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){ 
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }
  
  #readies for saving
  write.phy(aligned.set, file= paste0(align.dir, "/", trim.dir, "/", gsub(pattern = "\\..*", "", locus.names[i]), ".phy"), interleave = F)
  
}#end i loop
  
################################################################################################
################################################################################################
############## Step 6. Save Legacy loci separately                  ############################
################################################################################################
################################################################################################
################################################################################################

#Blast two probe set files together to get uce loci to keep
setwd(work.dir)

#Make blast database for the probe loci
system(paste("makeblastdb -in ", probe.file, " -parse_seqids -dbtype nucl ",
             " -out probe_blast_db", sep = ""))

#Matches samples to loci
system(paste("blastn -task dc-megablast -db probe_blast_db",
             " -query ", legacy.file, " -out leg_match.txt", 
             " -outfmt 6 -num_threads ", threads, sep = ""))

#headers for the blast db
headers<-c("qName", "tName", "pident", "matches", "misMatches", "gapopen", 
           "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

#Load in matches
match.data<-fread("leg_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
setnames(match.data, headers)
match.data<-match.data[match.data$evalue <= 0.05,]

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
probe.loci<-scanFa(FaFile(probe.file))
contigs<-scanFa(FaFile(legacy.file))
new.qsize<-width(contigs)[pmatch(filt.data$qName, names(contigs), duplicates.ok = T)]
filt.data[,qSize:=as.numeric(new.qsize)]
#Gets the sizes from the probes
new.tsize<-width(probe.loci)[pmatch(filt.data$tName, names(probe.loci), duplicates.ok = T)]
filt.data[,tSize:=as.numeric(new.tsize)]
#Removes matches that barely match to the full contig
filt.data<-filt.data[filt.data$matches >= filt.data$tSize *0.20,]

#Looks to remove duplicates
dup.match<-unique(filt.data[duplicated(filt.data$qName),]$qName)

keep.match<-c()
for (i in 1:length(dup.match)){
  
  sub.match<-filt.data[filt.data$qName %in% dup.match[i],]
  
  #Skips if they all match to each other
  if (length(unique(sub.match$qName)) == 1 && length(unique(sub.match$tName)) == 1){ next }
  
  #Looks for a pretty high bit score
  b.match<-sub.match[sub.match$bitscore/max(sub.match$bitscore) >= 0.6,]
  if (nrow(b.match) == 1){ 
    keep.match<-rbind(keep.match, b.match)
    next
  }
  
  c.match<-sub.match[sub.match$qSize/sub.match$tSize >= 0.6,]
  if (nrow(c.match) == 1){ 
    keep.match<-rbind(keep.match, c.match)
    next
  }
}

red.data<-filt.data[!filt.data$qName %in% dup.match,]
new.data<-rbind(red.data, keep.match)

#Gets the files to save separately
leg.names<-unique(new.data$tName)
leg.names<-paste(leg.names, ".phy", sep = "")
system("rm probe_blast* leg_match.txt")

#Sets up new directory for this stuff
trim.dir<-"legacy-markers_trimmed"
dir.create(paste0(align.dir, "/", trim.dir))
setwd(paste0(align.dir, "/exon-only_trimmed"))
file.names<-list.files(".")
locus.names<-leg.names[leg.names %in% file.names]

#Saves the new markers with new names 
for (i in 1:length(locus.names)){
  
  temp.data<-new.data[new.data$tName %in% gsub(".phy", "", locus.names[i]),]
  temp.data$qName<-gsub("_", "-", temp.data$qName)
  new.name<-paste0(temp.data$tName, "_", temp.data$qName)
  system(paste0("cp ", locus.names[i]," ", align.dir, "/", trim.dir, "/", new.name, ".phy"))

}

# END SCRIPT

