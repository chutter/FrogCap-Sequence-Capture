library(BioStrings)
library(PhyloCap)
library(data.table)
library(ggplot2)

#Locations
work.dir = "/Volumes/Armored/FrogCap_Anura_Seqcap"
assembly.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Final_Contigs/Anura_48"
target.file = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/01_Current-Version/FrogCap_Files/mtGenes.fa"
output.name = "mito-genomes"

#  settings
threads = 4
memory = 8
overwrite = TRUE
resume = FALSE
quiet = TRUE

# #program paths
blast.path = "/Users/chutter/miniconda3/bin"
bbmap.path = "/usr/local/bin"

####################################################################################
##### Do not edit below
setwd(work.dir)


#Add the slash character to path
if (is.null(blast.path) == FALSE){
  b.string = unlist(strsplit(blast.path, ""))
  if (b.string[length(b.string)] != "/") {
    blast.path = paste0(append(b.string, "/"), collapse = "")
  }#end if
} else { blast.path = "" }

#Add the slash character to path
if (is.null(bbmap.path) == FALSE){
  b.string = unlist(strsplit(bbmap.path, ""))
  if (b.string[length(b.string)] != "/") {
    bbmap.path = paste0(append(b.string, "/"), collapse = "")
  }#end if
} else { bbmap.path = "" }

#Initial checks
if (is.null(target.file) == T){ stop("A fasta file of targets to match to assembly contigs is needed.") }

#So I don't accidentally delete everything while testing resume
if (resume == TRUE & overwrite == TRUE){
  overwrite = FALSE
  stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
}

if (dir.exists(output.name) == TRUE) {
  if (overwrite == TRUE){
    system(paste0("rm -r ", output.name))
    dir.create(output.name)
  }
} else { dir.create(output.name) }

dir.create(paste0(output.name, "/sample-markers"))
dir.create(paste0(output.name, "/mitochondrial-markers"))
dir.create(paste0(output.name, "/species-raw-data"))

#Gets contig file names
file.names = list.files(assembly.directory)
ref.marker = readDNAStringSet(target.file)

#headers for the blast db
headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

#Matching and processing for each sample
all.seq = Biostrings::DNAStringSet()
for (i in 1:length(file.names)) {

  #Sets up working directories for each species
  sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
  species.dir = paste0(output.name, "/species-raw-data/", sample)

  #Creates species directory if none exists
  if (file.exists(species.dir) == F){ dir.create(species.dir) }

  #Checks if this has been done already
  if (resume == TRUE){
    if (file.exists(paste0(species.dir, "/", sample, "_", alignment.contig.name, "_matching-contigs.fa")) == T){ next }
  }#end

  #########################################################################
  #Part A: Blasting
  #########################################################################

  #Copies original contigs over
  system(paste0("cp ", assembly.directory, "/", file.names[i], " ",
                species.dir, "/", file.names[i]))

  # # DEDUPE almost exact duplicate removal
  system(paste0(bbmap.path, "dedupe.sh in=",assembly.directory, "/", file.names[i], " ordered=t overwrite=true ",
                " out=", species.dir, "/", sample, "_dedupe.fa", " minidentity=97"), ignore.stderr = quiet)

  #Make blast database for the probe loci
  system(paste0(blast.path, "makeblastdb -in ", species.dir, "/", sample, "_dedupe.fa",
                " -parse_seqids -dbtype nucl -out ", species.dir, "/", sample, "_nucl-blast_db"), ignore.stdout = quiet)

  #Removes extra file
  system(paste0("rm ", species.dir, "/", file.names[i]))

  #Matches samples to loci
  system(paste0(blast.path, "blastn -task dc-megablast -db ", species.dir, "/", sample, "_nucl-blast_db -evalue 0.001",
                " -query ", target.file, " -out ", species.dir, "/", sample, "_target-blast-match.txt",
                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                " -num_threads ", threads))

  #Need to load in transcriptome for each species and take the matching transcripts to the database
  system(paste0("rm ", species.dir, "/*nucl-blast_db*"))

  #Loads in match data
  match.data = data.table::fread(paste0(species.dir, "/", sample, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

  if (nrow(match.data) == 0) {
    print(paste0(sample, " had no matches. Skipping"))
    next }

  data.table::setnames(match.data, headers)
  filt.data = match.data[match.data$matches >= 20,]

  if (nrow(filt.data) == 0) {
    print(paste0(sample, " had no matches. Skipping"))
    next }

  #Sorting: exon name, contig name, bitscore higher first, evalue
  data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

  #Make sure the hit is greater than 50% of the reference length
  filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

  #Reads in contigs
  contigs = Biostrings::readDNAStringSet(paste0(species.dir, "/", sample, "_dedupe.fa"), format = "fasta")
  names(contigs) = gsub(" .*", "", names(contigs))

  # Goes through and cuts every locus
  marker.names = unique(filt.data$qName)

  sample.seq = Biostrings::DNAStringSet()
  for (j in 1:length(marker.names)){

    temp.data = filt.data[filt.data$qName == marker.names[j],]
    temp.data = temp.data[temp.data$bitscore == max(temp.data$bitscore),][1]
    temp.contig = contigs[names(contigs) %in% temp.data$tName]

    start = min(temp.data$tStart, temp.data$tEnd)
    end = max(temp.data$tStart, temp.data$tEnd)
    temp.seq = Biostrings::subseq(x = temp.contig, start = start, end = end)
    names(temp.seq) = paste0(sample, "_|_", marker.names[j] )
    sample.seq = append(sample.seq, temp.seq)

  }#end j loop

  #Save separate sample and separate marker files
  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(sample.seq))
  writeFasta(sequences = final.loci, names = names(final.loci),
             paste0(species.dir, "/", sample, "_mitochondrial-markers.fa"), nbchar = 1000000, as.string = T)

  all.seq = append(all.seq, sample.seq)

}#end i loop

#Saves as separate files
sample.names = unique(gsub("_\\|_.*", "", names(all.seq)))
marker.names = unique(gsub(".*_\\|_", "", names(all.seq)))

header.data = c("sample", "marker_name", "length", "ref_length", "completeness")
save.data = data.table::data.table(matrix(as.double(0), nrow = 0, ncol = length(header.data)))
data.table::setnames(save.data, header.data)
save.data[, sample:=as.character(sample)]

for (i in 1:length(marker.names)){

  temp.marker = all.seq[gsub(".*_\\|_", "", names(all.seq)) %in% marker.names[i]]
  temp.ref = ref.marker[names(ref.marker) %in% marker.names[i]]

  #header data
  header.data = c("sample", "marker_name", "length", "ref_length", "completeness")
  temp.data = data.table::data.table(matrix(as.double(0), nrow = length(temp.marker), ncol = length(header.data)))
  data.table::setnames(temp.data, header.data)
  temp.data[, sample:=as.character(sample)]

  #Gets the saved matching targets
  data.table::set(temp.data, j = match("sample", header.data), value = gsub("_\\|_.*", "", names(temp.marker)) )
  data.table::set(temp.data, j = match("marker_name", header.data), value = gsub(".*_\\|_", "", names(temp.marker)) )
  data.table::set(temp.data, j = match("length", header.data), value = width(temp.marker) )
  data.table::set(temp.data, j = match("ref_length", header.data), value = width(temp.ref) )
  data.table::set(temp.data, j = match("completeness", header.data), value =width(temp.marker) / width(temp.ref) )

  save.data = rbind(save.data, temp.data)

  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(temp.marker))
  writeFasta(sequences = final.loci, names = names(final.loci),
             paste0(output.name, "/mitochondrial-markers/", marker.names[i], ".fa"), nbchar = 1000000, as.string = T)

}#end i loop

#Finds probes that match to two or more contigs
final.loci = as.list(as.character(all.seq))
writeFasta(sequences = final.loci, names = names(final.loci),
           paste0(output.name, "/all-marker_unaligned.fa"), nbchar = 1000000, as.string = T)

#Saves combined, final dataset
save.data$completeness[save.data$completeness >= 1] = 1
write.csv(save.data, file = paste0(output.name, "/sample-marker_completeness-raw.csv"), row.names = F)


#### Make a heatmap
header.data = c("sample", marker.names)
heatmap.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
data.table::setnames(heatmap.data, header.data)
heatmap.data[, sample:=as.character(sample.names)]

for (i in 1:length(marker.names)){

  marker.data = save.data[save.data$marker_name %in% marker.names[i],]

  data.table::set(heatmap.data, i = match(marker.data$sample, sample.names),
                  j = as.integer(match(marker.names[i], marker.names)+1),
                  value = marker.data$completeness)
}


heatmap.plot = as.matrix(heatmap.data[,2:ncol(heatmap.data)])
rownames(heatmap.plot) = heatmap.data$sample
heatmap(heatmap.plot, Colv = NULL, Rowv = NULL, scale = "column", col = c("blue", "yellow", "orange", "red"))

###GGPLOT VERSION
plot.data = save.data[grep("tRNA-", save.data$marker_name, invert = T),]
ggplot(plot.data, aes(sample, marker_name, fill = completeness)) + geom_tile(aes(fill = completeness)) +
  scale_fill_viridis_b(option = "magma")

#Makes some plots
#temp.mean = unlist(lapply( split(save.data$completeness, save.data$sample), mean))
#temp.mean = unlist(lapply( split(save.data$completeness, save.data$marker_name), mean))

temp.mean = unlist(lapply( split(plot.data$completeness, plot.data$marker_name), mean))
mean(temp.mean)


