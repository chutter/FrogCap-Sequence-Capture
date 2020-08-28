#Required R packages
library(seqinr)
library(data.table)
library(Biostrings)

options(stringsAsFactors = FALSE)

###############################################################################
########   Step 1: Create locus summary dataset                #################
###############################################################################
###############################################################################
###############################################################################

###### Directory setup ########

work.dir = "/Full/Directory/Path/to/Custom-Probe-Set"
output.file ="Custom"

#The two metadata files should be placed in the working directory 
tax.data.file = "Taxonomy_dataset.csv"
marker.data.file = "Marker_Bait_dataset.csv"

#################################################
#Features you would like the probe set to have
#################################################

probe.size = 20000 #Number of probes, maximum is currently 20k for custom kits. use full kit
uce.include = TRUE #Other option FALSE, no UCEs
most.informative = TRUE #FALSE = randomly selects markers rather than by information content
min.exon.length = 200 #The minimum exon length
max.exon.length = 100000 #The maximum exon length

#Select a family or subfamily from the two lists below
clade.select = "Microhylidae"

#You can also select two or more taxonomic groups [overrides previous clade select]
#The program will select markers present in BOTH. 
#Must both be the same phylogenetic scale (Family / Subfamily)
#Cannot be from different Superfamilies (Ranoidea / Hyloidea)
clade.select = c("Mantellidae", "Rhacophoridae")

#Family group, options include:
#Ranoidea
#[1] "Mantellidae"         "Pyxicephalidae"      "Ranidae"            
#[4] "Microhylidae"        "Arthroleptidae"      "Ceratobatrachidae"  
#[7] "Hyperoliidae"        "Dicroglossidae"      "Rhacophoridae"      
#[10] "Odontobatrachidae"   

#Hyloidea
#"Aromobatidae"        "Bufonidae"          
#[13] "Craugastoridae"      "Salamander"          "Ceratophryidae"     
#[16] "Dendrobatidae"       "Eluetherodactylidae" "Hemiphractidae"     
#[19] "Centrolenidae"       "Hylidae"             "Megophryidae"       
#[22] "Leptodactylidae"     "Scaphiopodidae"      "Telmatobiidae" 


#Submfamily group, options include: 

#Ranoidea
#[1] "Laliostominae"       "Pyxicephalidae"      "Ranidae"            
#[4] "Cophylinae"          "Arthroleptinae"      "Boophinae"          
#[7] "Ceratobatrachinae"   "Gastrophryninae"     "Hyperoliinae"       
#[10] "Dicroglossinae"      "Kassininae"          "Astylosterninae"    
#[13] "Microhylinae"        "Rhacophorinae"       "Odontobatrachidae"  

#Hyloidea
#[16] "Asterophryinae"      "Scaphiophryninae"    "Aromobatidae"       
#[19] "Bufonidae"           "Holoadeninae"        "Salamander"         
#[22] "Ceratophryidae"      "Craugastorinae"      "Dendrobatidae"      
#[25] "Eluetherodactylidae" "Hemiphractidae"      "Hyalinobatrachinae" 
#[28] "Hylinae"             "Leptobrachiinae"     "Leptodactylinae"    
#[31] "Pelodryadinae"       "Lophyohylinae"       "Leiuperinae"        
#[34] "Ceuthomantinae"      "Scaphiopodidae"      "Telmatobiidae"      
#[37] "Centroleninae" 


#######################################################
#######################################################
#######################################################
#######################################################
#### Script running, do not edit below
#######################################################

setwd(work.dir)

#Reads in and modifies data
tax.data = fread(tax.data.file)
tax.data[,V1:=NULL]
marker.data = fread(marker.data.file)
marker.data$exon_length = nchar(marker.data$Marker_Seq)

#Reduces based on taxonomy selected
tax.red = tax.data[tax.data$Family %in% clade.select,]

if (nrow(tax.red) == 0){ tax.red = tax.data[tax.data$Subfamily %in% clade.select,] }
if (nrow(tax.red) == 0){stop("Problem with taxaonomic group selected.") }

### Pick markers where samples exist
clade.markers = colSums(tax.red[,6:ncol(tax.red)])
clade.markers = clade.markers[clade.markers >= 1]

if (length(clade.select) >= 2){
  keep.markers = c()
  for (i in 1:length(clade.select)){
    tax.temp = tax.red[tax.red$Family == clade.select[i],]
    if (nrow(tax.temp) ==0){ tax.temp = tax.red[tax.red$Subfamily == clade.select[i],] }
    temp.markers = colSums(tax.temp[,6:ncol(tax.temp)])
    temp.markers = temp.markers[temp.markers >= 1]
    keep.markers = append(temp.markers, keep.markers)
  }#end for loop
  
  clade.markers = table(names(keep.markers))
  clade.markers = clade.markers[clade.markers >= length(clade.select)]
}#end if statement

#Reduce down the  markers based on taxonomy
red.marker = marker.data[marker.data$ProbeSet %in% unique(tax.red$Probe_Set),]
red.marker = red.marker[red.marker$Marker_Name %in% names(clade.markers),]

if (sum(red.marker$BaitCount) <= probe.size){ stop("Not enough baits, consider adding more taxonomic groups.")}

#Reduce down the markers based on other characteristics

#Excludes UCEs
if (uce.include == FALSE){
  red.marker = red.marker[grep("uce", red.marker$Marker_Name, invert = T),]
}
if (sum(red.marker$BaitCount) <= probe.size){ stop("Not enough baits, consider modifying UCE inclusion.")}

#Selects the markers based on selected exon lengths
red.marker = red.marker[red.marker$exon_length >= min.exon.length,]
red.marker = red.marker[red.marker$exon_length <= max.exon.length,]
if (sum(red.marker$BaitCount) <= probe.size){ stop("Not enough baits, considering modifying exon length.")}

#Picks the most informative markers if selected
if (most.informative == TRUE){
  #Orders by informative sites
  order.marker = red.marker[order(red.marker$prop_pis, decreasing = T),]
  #While loop that adds markers until threshold probe set size is reached
  total.baits = 0
  custom.markers = c()
  count = 1
  while(total.baits <= probe.size){
    #Keeps adding markers until condition met
    custom.markers = rbind(custom.markers, order.marker[count,])
    #Checks counts
    count = count+1
    total.baits = sum(custom.markers$BaitCount)
  }
  
  #Finally, cuts out least informative marker that is over the probe threshold
  to.remove = total.baits - probe.size
  remove.markers = custom.markers[custom.markers$BaitCount == to.remove,]
  #Finalize
  custom.markers = custom.markers[custom.markers$Marker_Name != remove.markers$Marker_Name[nrow(remove.markers)],]
} #end most.informative if

#Picks markers randomly if informativeness does not matter
if (most.informative == FALSE){
  #Randomly orders dataset
  order.marker = red.marker[sample(nrow(red.marker)),]
  
  #While loop that adds markers until threshold probe set size is reached
  total.baits = 0
  custom.markers = c()
  count = 1
  while(total.baits <= probe.size){
    #Keeps adding markers until condition met
    custom.markers = rbind(custom.markers, order.marker[count,])
    #Checks counts
    count = count+1
    total.baits = sum(custom.markers$BaitCount)
  } #end while
  
  #Finally, randomly cut out markers over threshold to make even amount
  to.remove = total.baits - probe.size
  remove.markers = custom.markers[custom.markers$BaitCount == to.remove,]
  #Finalize
  custom.markers = custom.markers[custom.markers$Marker_Name != remove.markers$Marker_Name[nrow(remove.markers)],]
} #end most.informative if

#Save the collection of output files
######################################3

#Save the marker sequence
save.markers = DNAStringSet(as.character(custom.markers$Marker_Seq))
names(save.markers) = custom.markers$Marker_Name
final.markers = as.list(as.character(save.markers))
write.fasta(sequences = final.markers, names = names(final.markers), 
            paste0(output.file, "_Markers.fa"), nbchar = 1000000, as.string = T)

#Save the probe sequences, need to pull out of metadata
save.probes = DNAStringSet()
for (i in 1:nrow(custom.markers)){
  #Loops through each marker and saves probes separately
  temp.marker = custom.markers[i,]
  temp.probes = as.character(temp.marker$Bait_Seqs)
  split.probes = unlist(strsplit(temp.probes, split = "\\+"))
  new.probes = DNAStringSet(as.character(split.probes))
  names(new.probes) = paste0(temp.marker$Marker_Name, "_|_", rep(1:length(new.probes)))
  save.probes = append(save.probes, new.probes)
}  

#Saves them as fasta
final.probes = as.list(as.character(save.probes))
write.fasta(sequences = final.probes, names = names(final.probes), 
            paste0(output.file, "_Baits.fa"), nbchar = 1000000, as.string = T)

#Save the metadata
write.table(custom.markers, file = paste0(output.file, "_metadata.txt"), row.names = F)


####################################
### END SCRIPT
####################################
####################################
####################################
####################################
