library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

options(stringsAsFactors = FALSE)

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

###############################################################################
###############################################################################
######################Step 1 pull out loci and align  #########################
###############################################################################
###############################################################################

###### PARAMETER SETUP ####
threads<-"8"
min.taxa<-5 #min number to keep an alignment

#Join dirs
out.dir<-"Combine_Alignments"
work.dir<-"/Users/chutter/Dropbox/Research/WIP/Anura_Seqcap/Post_Processing"
align.hyloidea<-"/Users/chutter/Dropbox/Research/WIP/Anura_Seqcap/Post_Processing/Alignments/Alignments_hyloidea_broad"
align.ranoidea<-"/Users/chutter/Dropbox/Research/WIP/Anura_Seqcap/Post_Processing/Alignments/Alignments_ranoidea_broad"
probe.file<-"/Users/chutter/Dropbox/Research/WIP/FrogCap_Pipeline/Source_Files/Full_Ranoidea_Loci_May23.fa"

#Cluster dirs
#work.dir<-"/home/c111h652/scratch/Anura_Seqcap"
#probe.file<-"/home/c111h652/scratch/Anura_Seqcap/Full_Ranoidea_Loci_May23.fa"
#out.dir<-"/home/c111h652/scratch/Anura_Seqcap/Alignments_ranoidea_broad"
#species.loci<-"ranoidea_broad_24_contigs.fa"


#START
#########################
setwd(work.dir)
dir.create(out.dir)
setwd(paste0(work.dir, "/", out.dir))

#Gets locus names
dir.create("exon_only")
dir.create("exon_intron")

#Find alignments in common 
hyl.align<-list.files(paste0(align.hyloidea, "/all_untrimmed"))
ran.align<-list.files(paste0(align.ranoidea, "/all_untrimmed"))
share.align<-ran.align[ran.align %in% hyl.align]
bait.loci<-scanFa(FaFile(probe.file))  # loads up fasta file

#Loops through each locus and writes each species to end of file
for (i in 5079:length(share.align)){
  
  #STEP 1: Gather files to combine
  ##############
  hyl.align<-DNAStringSet(readAAMultipleAlignment(file = paste0(align.hyloidea,"/all_untrimmed/",share.align[i]), format = "phylip"))
  ran.align<-DNAStringSet(readAAMultipleAlignment(file = paste0(align.ranoidea,"/all_untrimmed/",share.align[i]), format = "phylip"))
  
  #Gets reference locus
  ref.locus<-bait.loci[names(bait.loci) %in% gsub(".phy$", "", share.align[i])]
  names(ref.locus)<-paste("Reference_Locus")

  #STEP 2: Runs MAFFT to add
  ##############
  temp.align<-run.mafft(unaligned.contigs = ran.align, add.contigs = hyl.align, rev.dir = T,
                       algorithm = "add", save.name = gsub(".phy$", "",share.align[i]), delete.files = T)
  
  alignment<-run.mafft(unaligned.contigs = temp.align, add.contigs = ref.locus, rev.dir = T,
                       algorithm = "add", save.name = gsub(".phy$", "",share.align[i]), delete.files = F)
  
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
    #Deletes old files
    system(paste("rm ", gsub(".phy$", "",share.align[i]), "_align.fa ", sep = ""))
    print(paste(share.align[i], " had too few taxa", sep = ""))
    next }
  
  ### realign if bad seqs removed
  if (length(bad.seqs) != 0){
    
    #Checks if its the entire aligment
    check.ran<-names(ran.align)[!names(ran.align) %in% bad.seqs]
    check.hyl<-names(hyl.align)[!names(hyl.align) %in% bad.seqs]
    if(length(check.hyl) <= 3 || length(check.ran <= 3)){ 
      system(paste("rm ", gsub(".phy$", "",share.align[i]), "_align.fa ", sep = ""))
      print(paste(share.align[i], " locus does not align betweetn ranoidea/hyloidea"))
      next
      }#end if
    
    #Runs mafft
    alignment<-run.mafft(rem.align, save.name = share.align[i], threads = threads, delete.files = F)  
    reversed<-names(alignment)[grep(pattern = "_R_", names(alignment))]
    if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment<-reverseComplement(alignment) }
    names(alignment)<-gsub(pattern = "_R_", replacement = "", x = names(alignment))
  } # end bad.seqs if
  
  ##############
  #STEP 4: Save UCE data as the same for both intron and exon + intron datasets since they don't have this
  ##############
  
  #Detects from the name
  if (length(grep(pattern = "uce", x=share.align[i])) == 1){
    #writes alignment
    red.align<-alignment[!names(alignment) %in% "Reference_Locus"]
    new.align<-strsplit(as.character(red.align), "")
    mat.align<-lapply(new.align, tolower)
    m.align<-as.matrix(as.DNAbin(mat.align))

    #readies for saving
    write.phy(m.align, file=paste0("exon_intron/", share.align[i]), interleave = F)
    write.phy(m.align, file=paste0("exon_only/", share.align[i]), interleave = F)
    
    #Deletes old files
    system(paste("rm ", gsub(".phy$", "",share.align[i]), "_align.fa ", sep = ""))
    print(paste(share.align[i], " UCE saved"))
    next
  }
  
  ##############
  #STEP 5: Make trimmed exon alignments
  ############## 
  
  #Removes the edge gaps
  ref.aligned<-as.character(alignment['Reference_Locus'])
  not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start<-min(not.gaps)
  ref.finish<-max(not.gaps)
  temp.intron<-subseq(alignment, ref.start, ref.finish)
  
  #Saves prelim exon file
  save.align<-temp.intron[!names(temp.intron) %in% "Reference_Locus"]
  new.align<-strsplit(as.character(save.align), "")
  mat.align<-lapply(new.align, tolower)
  aligned.set<-as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(aligned.set, file=paste0("exon_only/", share.align[i]), interleave = F)
  
  ##############
  #STEP 6: Make trimmed intron + exon alignments alignments
  ############## 
  
  #writes alignment
  red.align<-alignment[!names(alignment) %in% "Reference_Locus"]
  new.align<-strsplit(as.character(red.align), "")
  mat.align<-lapply(new.align, tolower)
  m.align<-as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(m.align, file=paste0("exon_intron/", share.align[i], sep = ""), interleave = F)
  
  #Deletes old files
  system(paste0("rm ", gsub(".phy$", "",share.align[i]), "_align.fa "))
  print(paste0(gsub(".phy$", "",share.align[i]), " Locus saved"))

}# end big i loop  


#END SCRIPT

