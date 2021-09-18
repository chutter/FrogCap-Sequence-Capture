#' @title summary.sampleSpecificity
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

summary.sampleSpecificity = function(read.directory = NULL,
                                     sub.directory = NULL,
                                     target.fasta = NULL,
                                     output.name = "sample-specificity",
                                     sample.groups = NULL,
                                     threads = 1,
                                     memory = 1,
                                     overwrite = FALSE,
                                     quiet = FALSE,
                                     bwa.path = NULL,
                                     picard.path = NULL,
                                     samtools.path = NULL) {


  # #Directories and stuff
  #  ran.in<-Biostrings::readDNAStringSet("/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Final_Ranoidea_Probe-Loci_Sept28.fa")
  #  red.in<-Biostrings::readDNAStringSet("/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Final_RedRan_Probe-Loci_Sept28.fa")
  #  hyl.in<-Biostrings::readDNAStringSet("/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Final_Hyloidea_Probe-Loci_Sept28.fa")
  #
  #  ran.in = ran.in[duplicated(names(ran.in)) == FALSE]
  #  red.in = red.in[duplicated(names(red.in)) == FALSE]
  #  hyl.in = hyl.in[duplicated(names(hyl.in)) == FALSE]
  #
  #  save.rownames = names(ran.in)
  #  write.align = as.list(as.character(ran.in))
  #
  #  #Creates random name and saves it
  #  writeFasta(sequences = write.align,
  #             names = names(write.align),
  #             file.out = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/RedRan_Aug18-2021.fa",
  #             nbchar = 1000000,
  #             as.string = T)
  #
  #
  #  save.rownames = names(red.in)
  #  write.align = as.list(as.character(red.in))
  #
  #  #Creates random name and saves it
  #  writeFasta(sequences = write.align,
  #             names = names(write.align),
  #             file.out = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Ranoidea_Aug18-2021.fa",
  #             nbchar = 1000000,
  #             as.string = T)
  #
  #
  #  save.rownames = names(hyl.in)
  #  write.align = as.list(as.character(hyl.in))
  #
  #  #Creates random name and saves it
  #  writeFasta(sequences = write.align,
  #             names = names(write.align),
  #             file.out = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Hyloidea_Aug18-2021.fa",
  #             nbchar = 1000000,
  #             as.string = T)
  #
  #

  # ran.loci = ran.in[names(ran.in) %in% names(hyl.in)]
  # hyl.loci = hyl.in[names(hyl.in) %in% names(ran.in)]
  # red.loci = red.in[names(red.in) %in% names(ran.in)]
  #
  #
  # #Finds probes that match to two or more contigs
  # ran.loci = ran.loci[duplicated(names(ran.loci)) != T]
  # save.rownames = names(ran.loci)
  # write.align = as.list(as.character(ran.loci))
  #
  # #Creates random name and saves it
  # writeFasta(sequences = write.align,
  #            names = names(write.align),
  #            file.out = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Overlapping_Loci_Ran-Hyl_Aug18-2021.fa",
  #            nbchar = 1000000,
  #            as.string = T)

#
#   # #Debug
#   library(foreach)
#   setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
#   read.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Processed_Samples"
#   target.fasta = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Overlapping_Loci_Ran-Hyl_Aug18-2021.fa"
#   sub.directory = "cleaned-reads-snp"
#   sample.groups = "sample_probeset_group.csv"
#   output.name = "Analyses/sample-specificity"
#   samtools.path = "/Users/chutter/miniconda3/bin"
#   bwa.path = "/usr/local/bin"
#   picard.path = "/Users/chutter/miniconda3/bin"
#   bbmap.path = "/usr/local/bin"
#   overwrite = TRUE
#
#
  read.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Processed_Samples"
  sub.directory = "cleaned-reads-snp"
  target.shared = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Overlapping_Loci_Ran-Hyl_Aug18-2021.fa"
  target.ranoidea = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Ranoidea_Aug18-2021.fa"
  target.redran = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/RedRan_Aug18-2021.fa"
  target.hyloidea = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Hyloidea_Aug18-2021.fa"
  sample.groups = read.csv("sample_probeset_group.csv", header = T)
  ranoidea.groups = sample.groups[sample.groups[,1] %in% "Ranoidea-V1",]
  hyloidea.groups = sample.groups[sample.groups[,1] %in% "Hyloidea-V1",]
  redran.groups = sample.groups[sample.groups[,1] %in% "RedRan",]

  samtools.path = "/Users/chutter/miniconda3/bin"
  bwa.path = "/usr/local/bin"
  picard.path = "/Users/chutter/miniconda3/bin"
  bbmap.path = "/usr/local/bin"

  output.name = "Analyses/sample-specificity"

  target.fasta = target.redran
  sample.groups = redran.groups
  quiet = TRUE
  threads = 6
  memory = 6


  ##### Program path check
  ####################################################################
  #Same adds to bbmap path
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Same adds to bbmap path
  if (is.null(picard.path) == FALSE){
    b.string = unlist(strsplit(picard.path, ""))
    if (b.string[length(b.string)] != "/") {
      picard.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { picard.path = "" }
  ####################################################################

  #Quick checks
  if (is.null(read.directory) == TRUE){ stop("Please provide input directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.name) == F){ dir.create(output.name) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.name))
      dir.create(output.name)
    }
  }#end else

  #  #Read in sample data **** sample is run twice?!
  reads = list.files(read.directory, recursive = T, full.names = T)
  if (is.null(sub.directory) != TRUE) {
    reads = reads[grep(paste0(sub.directory, "/"), reads)]
    sample.names = gsub(paste0("/", sub.directory, "/.*"), "", reads)
    sample.names = unique(gsub(paste0(read.directory, "/"), "", sample.names))
  } else {
    sample.names = list.files(read.directory, recursive = F, full.names = F)
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #Read in group data and the target markers
  if (class(sample.groups) == "character"){
    group.data = read.csv(sample.groups, header = T)
    } else {
    group.data = sample.groups
    }#end else

  sample.names = sample.names[sample.names %in% group.data[,2]]
  target.markers = Biostrings::readDNAStringSet(target.fasta)

  #Specificity refers to the percentage of cleaned reads that can be mapped back to the target markers.

  #Sets up data to collect
  header.all = c("dataset", "sample", "mapped_reads", "unmapped_mates",  "total_unmapped_reads",
                  "total_read_pairs",  "median_marker_rpkm", "mean_marker_rpkm", "sample_specificity")
  all.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.all)))
  data.table::setnames(all.data, header.all)
  all.data[, sample:=as.character(sample)]
  all.data[, dataset:=as.character(dataset)]

  #Loops through each locus and does operations on them
  for (i in 1:length(sample.names)){

    #Reads in the contig file to check
    dir.create(paste0(output.name, "/", sample.names[i]))
    #     sample.contigs = paste0(read.directory, "/", sample.names[i], "/assembled-contigs/", sample.names[i], "_orthologs.fa")

    #Sets up data to collect
    header.data = c("dataset","sample", "marker", "mapped_reads", "unmapped_mates",
                    "marker_total_reads", "total_unmapped", "marker_rpkm")
    collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(target.markers), ncol = length(header.data)))
    data.table::setnames(collect.data, header.data)
    collect.data[, sample:=as.character(sample)]
    collect.data[, dataset:=as.character(dataset)]
    collect.data[, marker:=as.character(marker)]

    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }

    sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    for (j in 1:length(sample.reads)){

      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      #Checks the Sample column in case already renamed
      if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
      #Returns an error if reads are not found
      if (length(lane.reads) == 0 ){
        stop(sample.reads[j], " does not have any reads present for files ")
      } #end if statement

      #Copies ortholog file to act as a reference
      #       system(paste0("cp ", sample.contigs, " ", output.name, "/", sample.names[i], "/sample-reference.fa"))

      system(paste0("cp ", target.fasta, " ", output.name, "/", sample.names[i], "/sample-reference.fa"))
      system(paste0(bwa.path, "bwa index ", output.name, "/", sample.names[i], "/sample-reference.fa"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #Creates a bam alignment file of reads mapped to reference
      system(paste0(bwa.path, "bwa mem -M -t ", threads, " ",
                    output.name, "/", sample.names[i], "/sample-reference.fa ",
                    lane.reads[1], " ", lane.reads[2],
                    " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
                    " -o ", output.name, "/", sample.names[i], "/paired.bam  -"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      system(paste0(picard.path, "picard -Xmx", memory, "g",
                    " MarkDuplicates INPUT=", output.name, "/", sample.names[i], "/paired.bam",
                    " OUTPUT=", output.name, "/", sample.names[i], "/dedup_paired.bam",
                    " METRICS_FILE=", output.name, "/", sample.names[i], "/metrics.txt",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      system(paste0("mv ", output.name, "/", sample.names[i], "/metrics.txt ",
                    output.name, "/", sample.names[i], "/duplication_stats.txt"))

      #HERe again
      system(paste0(picard.path, "picard -Xmx", memory, "g",
                    " BuildBamIndex INPUT=", output.name, "/", sample.names[i], "/paired.bam",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #Gets the genome coverage using bedtools
      system(paste0(samtools.path, "samtools idxstats --threads ", threads, " ",
                    output.name, "/", sample.names[i], "/paired.bam > ",
                    output.name, "/", sample.names[i], "/samtools_idxstats.txt"))

      #system(paste0("rm ", output.name, "/", sample.names[i], "/paired.bam"))

      #IDX data
      i.headers = c("marker", "length", "mapped_reads", "unmapped_mates")
      idx.stats = data.table::fread(file = paste0(output.name, "/", sample.names[i], "/samtools_idxstats.txt"))
      data.table::setnames(idx.stats, i.headers)
      idx.stats = idx.stats[idx.stats$marker != "*",]

      #Gets the number of reads that mapped to the reference, also total reads
      #Total number of un-mapped reads
      mapped.all = sum(idx.stats$mapped_reads)
      mapped.all = as.numeric(system(paste0("samtools view -c -F 4 ",
                                             output.name, "/", sample.names[i], "/paired.bam"), intern = T))
      unmapped.all = as.numeric(system(paste0("samtools view -c -f 4 ",
                                              output.name, "/", sample.names[i], "/paired.bam"), intern = T))

      #Goes through each locus and calculates stats
      locus.names = unique(idx.stats$marker)
      group.name = group.data[group.data[,2] %in% sample.names[i],][,1]
      locus.rpkm = idx.stats$mapped_reads / (idx.stats$length/1000 * mapped.all/1000000)

      #Starter stats
      data.table::set(collect.data, j = match("dataset", header.data), value = group.name )
      data.table::set(collect.data, j = match("sample", header.data), value = sample.names[i] )
      data.table::set(collect.data, j = match("marker", header.data), value = locus.names )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("mapped_reads", header.data), value = idx.stats$mapped_reads )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("unmapped_mates", header.data), value = idx.stats$unmapped_mates )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("marker_total_reads", header.data), value = (idx.stats$unmapped_mates+idx.stats$mapped_reads) )
      data.table::set(collect.data, j = match("total_unmapped", header.data), value = unmapped.all )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("marker_rpkm", header.data), value = locus.rpkm )

    } #end j loop

    write.table(collect.data, file = paste0(output.name, "/", sample.names[i], "/sample_specificity_raw-data.txt"), sep = "\t", row.names = F)

    #Collects the overall summary data
    collect.data = collect.data[collect.data$mapped_reads != 0,]
    data.table::set(all.data, i = as.integer(i), j = match("dataset", header.all), value = group.name )
    data.table::set(all.data, i = as.integer(i), j = match("sample", header.all), value = sample.names[i] )
    data.table::set(all.data, i = as.integer(i), j = match("mapped_reads", header.all), value = sum(collect.data$mapped_reads) )
    data.table::set(all.data, i = as.integer(i), j = match("unmapped_mates", header.all), value = sum(collect.data$unmapped_mates) )
    data.table::set(all.data, i = as.integer(i), j = match("total_unmapped_reads", header.all), value = unmapped.all )
    data.table::set(all.data, i = as.integer(i), j = match("total_read_pairs", header.all), value = sum(collect.data$marker_total_reads)+unmapped.all )
    data.table::set(all.data, i = as.integer(i), j = match("mean_marker_rpkm", header.all), value = mean(collect.data$marker_rpkm) )
    data.table::set(all.data, i = as.integer(i), j = match("median_marker_rpkm", header.all), value = median(collect.data$marker_rpkm) )
    data.table::set(all.data, i = as.integer(i), j = match("sample_specificity", header.all), value = sum(collect.data$marker_total_reads)/(sum(collect.data$marker_total_reads)+unmapped.all) )

    print(paste0(sample.names[i], " Completed!"))

  }# end i loop

  write.table(all.data, file = paste0(output.name, "/sample-specificity_summary.txt"), sep = "\t", row.names = F)

  #### Print some textual summary here

}# end function

### END SCRIPT
