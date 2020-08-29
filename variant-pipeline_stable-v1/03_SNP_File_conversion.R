library(ape)
library(stringr)
library(data.table)
library(ShortRead)
library(Rsamtools)

#SNP packages
library(vcfR)
library(VariantAnnotation)
library(TVTB)
library(snpStats)
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

#Filtering is not done through hard thresholds, but rather quality ranks
#Its all relative though, so if the data is bad, the SNPs will be bad. 
#Can use hard filtering to remove loci that are low quality 

###### parameter setup ########
threads = 2 # number of cores = uses 1 right now
mem.val = 8 #in GB

#Setup directories
out.dir<-"/your/work/directory/SNP-Analysis"
ref.dir = "/your/work/directory/Reference"

#From previous script
work.dir<-"/your/work/directory"
proc.dir<-"/your/work/directory/Processed_Samples"
align.dir<-"/your/work/directory/Alignments"
run.name<-"Output-Name"

#Setup directories
sample.file = "/your/work/directory/pop-1.csv"
work.dir = "/your/work/directory"
db.dir =  "/your/work/directory/SNP_Database"
dataset.name = "Name-your-dataset"

#remove.taxa = c() # If you don't want to remove any
remove.taxa = c("remove-species1", "remove-species2")

#Completed: 1
#Format options: 1: Structure, 2: G-PHOCS; 3: SVDQuartets, 4: Haplotype contigs
end.format = "Structure"

#Filters
filt.miss = 0.1 #proportion of missing data acceptable on a sample / site basis
filt.biall = FALSE #Use only sites that are biallelic ie. binary coding for certain programs
filt.poly = TRUE #Use only polymorphic sites? (yes recommended)
most.var = TRUE #Picks the most variable site from a locus, not sure if biases
resume = FALSE  #Resumes if crashes, bugs in VCF R memory allocation

#To do
#Add in population number
#Choose SNP that differentiates your populations the best 
#Remove taxa

############################################################################################
########### Part 2 #########################################################################
##### Converts from multisample VCF to Structure format  
############################################################################################

#Get loci names I guess
#probe.loci<-scanFa(FaFile(paste0(ref.dir, "/reference.fa")))
#loci.names<-names(probe.loci)

#Get same vcf files together
setwd(work.dir)

#Get multifile databases together
sample.data = read.csv(sample.file)
all.files = list.files(db.dir, recursive = F)
vcf.files = all.files[grep(".vcf$", all.files)]

#The second is structure format. I've attached the manual. See pages 4-7. Most of the extra data columns are not necessary. The only data needed is:
#1) Row 1: Marker names (this can just be a sequence of numbers from 1 - nloci)
#2) Row 2 onwards: 2 rows for each sample (diploid)
#3) After sample names, no other additional data is necessary so you can go straight into the genotyping (details in the manual; missing data is usually encoded with -9). 
#I've attached an example of a structure file from my amolops paper. That one actually has population designation in the second column (listed as 1, 2, 3) 
#before the genotypes (column 3 onwards). It would be ideal but not entirely necessary to write population designation into the script. It can be done post-hoc. 

#Sets up multiprocessing
#cl = makeCluster(threads)
#registerDoParallel(cl)
#mem.val<-floor(mem.val/threads)
#Loops through each locus and does operations on them
#foreach(i=1:length(loci.names), .packages = c("foreach")) %dopar% {

if (resume == TRUE){
  #Checks if file exists
  if(file.exists(paste0(dataset.name, "_structure-snps.str")) == T){
    final.locus = read.table(file = paste0(dataset.name, "_structure-snps.str"))
    colnames(final.locus) = final.locus[1,]
    final.locus = final.locus[-1,]
    start.pos = which(gsub(".vcf$", "", vcf.files) == colnames(final.locus)[ncol(final.locus)]) + 1
  #If not matching
  } else {
    start.pos = 1
    final.locus = data.frame(Label = sample.data$Sample, Pop = sample.data$Population) 
    final.locus = rbind(final.locus, final.locus)
    final.locus = final.locus[order(final.locus$Label),]
    if (length(remove.taxa) != 0){ final.locus = final.locus[!final.locus$Label %in% remove.taxa,] }
  }#end else
    
} else {
  #Sets up data
  start.pos = 1
  final.locus = data.frame(Label = sample.data$Sample, Pop = sample.data$Population) 
  final.locus = rbind(final.locus, final.locus)
  final.locus = final.locus[order(final.locus$Label),]
  if (length(remove.taxa) != 0){ final.locus = final.locus[!final.locus$Label %in% remove.taxa,] }
} #end 
  
  
#Loops through each locus and does operations on them
for (i in start.pos:length(vcf.files)){
  
  #Load in VCF file
 # vcf.locus = VariantAnnotation::readVcf(file = paste0(db.dir, "/", vcf.files[i]))
  #Gets genotypes
  #locus.gt = genotypeToSnpMatrix(vcf.locus, uncertain = F)
  
  #Maybe need
  #dna <- ape::read.dna(dna_file, format = "fasta")
  
  #processes and isolates high quality SNPs
  #vcf.locus = read.vcfR(paste0(db.dir, "/", vcf.files[i]), verbose = FALSE)
  
  
  vcf.locus = tryCatch(expr = {read.vcfR(paste0(db.dir, "/", vcf.files[i]), verbose = FALSE) },
                       error = function(x) { x<-data.frame(); return(x) })
  
  #Catches exceptions or bad files
  if (nrow(vcf.locus) == 0){ 
    print(paste0(vcf.files[i], " Could not be read. Something corrupted."))
    next 
    }
  
  #Also catches other types of bad files that load in
  if(nrow(vcf.locus@fix) == 0){ 
    print(paste0(vcf.files[i], " Had no initial variants. Skip."))
    next 
  }
  
  #Gts the data
  chrom = create.chromR(name='Supercontig', vcf=vcf.locus)
  filt.snps = masker(chrom, min_QUAL = 1)
  proc.snps = proc.chromR(filt.snps, verbose=F)

  bad.pos = which(is.na(proc.snps@var.info$MQ) == T)
  if (length(bad.pos) != 0){ proc.snps@var.info$mask[bad.pos] = FALSE }
  
  #Keeps sites that are heterozygous
  gt = extract.gt(vcf.locus)
  hets = is_het(gt)
  not.het = which(apply(hets, 1, function(x) sum(x == T))/ncol(hets) == 0)
  if (length(not.het) != 0){ proc.snps@var.info$mask[not.het] = FALSE }
  
  #Uses only biallelic sites
  if (filt.biall == T){
    #Gets the sites and masks them
    not.bi = which(is.biallelic(vcf.locus) == T)
    if (length(not.bi) != 0){ proc.snps@var.info$mask[not.bi] = FALSE }
  }#end if statement
  
  #Uses only polymorphic sites
  if (filt.poly == T){
    #Gets the sites and masks them
    not.po = which(is.polymorphic(vcf.locus) == F)
    if (length(not.po) != 0){ proc.snps@var.info$mask[not.po] = FALSE }
  }#end if statement
  
  #Filter out missing data
  md = extract.gt(proc.snps, element="GT", as.numeric=TRUE)
  not.good = which(apply(md, 1, function(x) sum(is.na(x)))/ncol(md) > filt.miss)
  if (length(not.good) != 0){ proc.snps@var.info$mask[not.good] = FALSE }
  
  #Ranks SNPs so the best can be chosen 
  gq = extract.gt(proc.snps, element="GQ", as.numeric=TRUE)
  dp = extract.gt(proc.snps, element="DP", as.numeric=TRUE)
  
  #Visualizes
  #chromoqc(proc.snps, dp.alpha=20)
  
 # par( mar = c(8,4,4,2) )
 # boxplot(gq, las=2, col=2:5, main="Genotype Quality (GQ)")
  
 # dp2 <- dp
  #dp2[ dp2 == 0 ] <- NA
 # boxplot(dp2, las=2, col=2:5, main="Sequence Depth (DP)", log="y")
 # abline(h=10^c(0:4), lty=3, col="#808080")
  
  #Finds optimal depth
  mids = apply(dp, MARGIN=2, median, na.rm=TRUE)
  dp2 = sweep(dp, MARGIN=2, mids, FUN="-")
  dp2 = abs(dp2)
  dp2 = -1 * dp2
  
  #Finds good genotypes
  gq2 = gq/100
  amins = abs(apply(dp2, MARGIN=2, min, na.rm = TRUE))
  dp2 = sweep(dp2, MARGIN=2, STATS = amins, FUN="+")
  dp2 = sweep(dp2, MARGIN=2, STATS = amins, FUN="/")
  #range(dp2, na.rm=TRUE)
  
  #par( mar = c(8,4,4,2) )
  #boxplot(dp2, las=2, col=2:5, main="Sequence Depth (DP)")
  
  #Gets a new score
  scores = dp2 + gq2
  scores = rowSums(scores, na.rm = TRUE)
  
  #hist(scores, col=4)
  
  #Rank variants with scores
  rank.snps = rank.variants.chromR(proc.snps, scores)
  #head(rank.snps@var.info)
  #rank.snps@var.info$mask[rank.snps@var.info$rank > 5] = FALSE
  
  #Select only heterozygous sites above a threshold
  #rank.snps@var.info$mask[rank.snps@var.info$He < filt.het] = FALSE
  het.sel = rank.snps@var.info[rank.snps@var.info$mask == T,]
  
  #If everythign is masked, no point to continue
  if (nrow(het.sel) == 0){
    print(paste0(vcf.files[i], " Had no good SNPs. Skipped."))
    next
  }
  
  #Filters to best 25% of the SNPS
  if (nrow(het.sel) >= 10){
    max.rank = (max(het.sel$rank) - min(het.sel$rank))/4
    rank.snps@var.info$mask[rank.snps@var.info$rank > max.rank] = FALSE
  } 
    
  #Extract best random genotype
  best.gt = extract.gt(rank.snps, element="GT", mask = T, as.numeric = F)

  #if there is one best.gt, it changes from matrix to character class.
  if (class(best.gt) == "character"){
    #Gets metadata back
    met.dat = rank.snps@var.info[rank.snps@var.info$mask == T,]
    temp.gt = data.frame(best.gt)
    colnames(temp.gt) = paste0(met.dat$CHROM, "_", met.dat$POS)
    best.gt = t(temp.gt)
  }#end outer if
    
  #Remove taxa
  if (length(remove.taxa) != 0){ best.gt = best.gt[,!colnames(best.gt) %in% remove.taxa] }
  
  #if there is one best.gt, it changes from matrix to character class.
  if (class(best.gt) == "character"){
    #Gets metadata back
    met.dat = rank.snps@var.info[rank.snps@var.info$mask == T,]
    temp.gt = data.frame(best.gt)
    colnames(temp.gt) = paste0(met.dat$CHROM, "_", met.dat$POS)
    best.gt = t(temp.gt)
  }#end outer if
  
  #Can fix if there are no genotypes, change threshold
  if(nrow(best.gt) == 0){ stop("Something bad happened. Let Carl know.") }
  
  #Save the data here in some sort of convertable format for later. 
  
  #Randomly pick the best site, maybe reduce down to most variable with a certain 
  site.het = apply(best.gt, MARGIN = 1, FUN = function (x) length(x[x!="0/0"])) / ncol(best.gt)
  sort.het = site.het[order(site.het, decreasing = T)]
  
  #Picks the most variable if you want
  if (most.var == T){ 
    var.site = sort.het[1]
    final.site = best.gt[rownames(best.gt) == names(var.site),]
  } #end most var
  
  if (most.var == F){
    #Sets up stuff
    stop = F
    count = 0
    while (stop == F) {
    
      #Randomly picks from the top 10% of variable SNPS
      red.het = sort.het[1:ceiling(length(sort.het)/10)]
      rand.site = red.het[sample(1:length(red.het), 1)]
      
      #Text break apart and create a data structure for the structure program 
      final.site = best.gt[rownames(best.gt) == names(rand.site),]
      
      if (length(unique(final.site)) <= 2){
        #break
        
        #Randomly picks from the top 10% of variable SNPS
        red.het = sort.het[1:ceiling(length(sort.het)/20)]
        rand.site = red.het[sample(1:length(red.het), 1)]
        
        #Text break apart and create a data structure for the structure program 
        final.site = best.gt[rownames(best.gt) == names(rand.site),]
        
        print("rerolling sites")
        
      } else {
        stop = T
      } #end if
      
      #stops if we get stuck in infinite loop
      count = count + 1
      if (count >= length(red.het)){ stop = T }#end if for count
      
    }#end while
    
    if (count >= length(red.het)){ 
      print(paste0(vcf.files[i], " Had no good SNPs because of taxa removal. Skipped."))
      next
    }#end if 
    
  } #end if for most.var = F
  

  #STRUCTURE FORMAT
  #1) Row 1: Marker names (this can just be a sequence of numbers from 1 - nloci)
  #2) Row 2 onwards: 2 rows for each sample (diploid)
  #3) After sample names, no other additional data is necessary so you can go straight into the genotyping (details in the manual; missing data is usually encoded with -9). 
  #I've attached an example of a structure file from my amolops paper. That one actually has population designation in the second column (listed as 1, 2, 3) 
  #before the genotypes (column 3 onwards). It would be ideal but not entirely necessary to write population designation into the script. It can be done post-hoc. 
  
  #Data frame structure
  #Label  Locus1 Locus2
  #Samp1    0     0
  #Samp1    1     0
  #Samp2    0     1
  #Samp2    0     1
  
  #Code missing data for structure
  final.site[which(is.na(final.site) == T)] = "9/9"
  
  #Split up the text codings into numeric row format
  temp.frame = unlist(strsplit(final.site, split = ""))
  temp.frame = temp.frame[temp.frame != "/"]
  temp.frame = temp.frame[temp.frame != "|"]
  names(temp.frame) = gsub(".$", "", names(temp.frame))
  locus.data = data.frame(Label = names(temp.frame), as.numeric(temp.frame), row.names = NULL)
  colnames(locus.data)[2] = gsub(".vcf$", "", vcf.files[i])
  
  #Missing data is -9 for structure 
  locus.data[locus.data == 9] = -9
  locus.data = locus.data[order(locus.data$Label),]
  
  #Adds to combined data
  final.locus = cbind(final.locus, locus.data[,2])
  names(final.locus)[ncol(final.locus)] = colnames(locus.data)[2]
  
  #Saves
  write.table(final.locus, paste0(dataset.name, "_structure-snps.str"), 
              row.names = F, sep = " ", quote = F, col.names = T)
  
}#end i loop  

colnames(final.locus)[3:length(colnames(final.locus))] = seq(1:(length(colnames(final.locus))-2))
final.locus$Label = gsub(".*_", "", final.locus$Label)

write.table(final.locus, paste0(dataset.name, "_structure-snps.str"), 
            row.names = F, sep = " ", quote = F, col.names = F)


  

#########################
###### END SCRIPT
#########################


  