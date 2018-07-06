library(ape)
library(seqinr)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)

options(stringsAsFactors = FALSE)

# Adapted di2multi function from the ape package to plot polytomies
# based on numeric node support values
# (di2multi does this based on edge lengths)
# Needs adjustment for unrooted trees as currently skips the first edge
di2multi4node <- function (phy, tol = 0.5) {
  if (is.null(phy$edge.length)) 
    stop("the tree has no branch length")
  if (is.na(as.numeric(phy$node.label[2])))
    stop("node labels can't be converted to numeric values")
  if (is.null(phy$node.label))
    stop("the tree has no node labels")
  ind <- which(phy$edge[, 2] > length(phy$tip.label))[as.numeric(phy$node.label[2:length(phy$node.label)]) < tol]
  n <- length(ind)
  if (!n) 
    return(phy)
  foo <- function(ancestor, des2del) {
    wh <- which(phy$edge[, 1] == des2del)
    for (k in wh) {
      if (phy$edge[k, 2] %in% node2del) 
        foo(ancestor, phy$edge[k, 2])
      else phy$edge[k, 1] <<- ancestor
    }
  }
  node2del <- phy$edge[ind, 2]
  anc <- phy$edge[ind, 1]
  for (i in 1:n) {
    if (anc[i] %in% node2del) 
      next
    foo(anc[i], node2del[i])
  }
  phy$edge <- phy$edge[-ind, ]
  phy$edge.length <- phy$edge.length[-ind]
  phy$Nnode <- phy$Nnode - n
  sel <- phy$edge > min(node2del)
  for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node2del < 
                                                           phy$edge[i])
  if (!is.null(phy$node.label)) 
    phy$node.label <- phy$node.label[-(node2del - length(phy$tip.label))]
  phy
}

###############################################################################
###############################################################################
######################Step 1 pull out loci and align  #########################
###############################################################################
###############################################################################

###### PARAMETER SETUP ####
threads<-"6"
min.taxa<-3 #min number to keep an alignment
output.file<-"microhylidae_trans_pf"

#Trees
work.dir<-"/Users/chutter/Dropbox/Research/WIP/Microhylidae_SeqCap/Tree_Analysis/Transcriptome/Astral"
tree.dir<-"/Users/chutter/Dropbox/Research/WIP/Microhylidae_SeqCap/Tree_Analysis/Transcriptome/Gene_Trees/IQTrees_trans_m"


##########################

setwd(tree.dir)
tree.files<-list.files(".")

#Loop through trees and collapse poorly supported nodes into polytomies
tree.list<-c()
for (i in 1:length(tree.files)){
  #read in tree file
  temp.tree<-read.tree(tree.files[i])
  if(length(temp.tree$node.label) == 0){
    write.tree(temp.tree, file = paste0(work.dir, "/", output.file, "_genetrees.tre"), append = T)
  } else {
    #replace empty with 0s
    temp.tree$node.label[temp.tree$node.label == ""]<-"0"
    #Find and collapse nodes with bl close to 0 from above 
    new.tree<-di2multi4node(temp.tree, tol = 10) 
    write.tree(new.tree, file = paste0(work.dir, "/", output.file, "_genetrees.tre"), append = T)
    
     #Add to list of trees
    #tree.list[[i]]<-new.tree
  }#end else
}#end if

#Save final tree file into astral folder
setwd(work.dir)
#Run Astral 
system(paste0("java -Xmx8g -jar /usr/local/Astral/astral.5.6.2.jar -i ", 
              output.file, "_genetrees.tre > ", 
              output.file, "_astral.tre"))








