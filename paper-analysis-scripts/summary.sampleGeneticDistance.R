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

summary.sampleGeneticDistance = function(alignment.directory = NULL,
                                         output.name = "sample-genetic-distance",
                                         target.type = c("reference", "consensus"),
                                         target.fasta = NULL,
                                         threads = 1,
                                         memory = 1,
                                         mafft.path = NULL,
                                         overwrite = FALSE,
                                         quiet = FALSE) {

  #Debug
  # setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
  #
  # ###*** Use the probes themselves!!
  #
  # #Alignment directoires
  # align.ranoidea = "/Users/chutter/Dropbox/Research/1_Main-Projects/2_Finished-Submitted/Hutter_etal_Anura_FrogCap/Post_Processing_New/Final_Alignments/Alignments_Ranoidea/all-markers_untrimmed"
  # align.hyloidea = "/Users/chutter/Dropbox/Research/1_Main-Projects/2_Finished-Submitted/Hutter_etal_Anura_FrogCap/Post_Processing_New/Final_Alignments/Alignments_Hyloidea/all-markers_untrimmed"
  # align.reduced = "/Users/chutter/Dropbox/Research/1_Main-Projects/2_Finished-Submitted/Hutter_etal_Anura_FrogCap/Post_Processing_New/Final_Alignments/Alignments_Reduced/all-markers_untrimmed"
  #
  # #Probe file directories
  # target.shared = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Overlapping_Loci_Ran-Hyl_Aug18-2021.fa"
  # target.ranoidea = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Ranoidea_Aug18-2021.fa"
  # target.reduced = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/RedRan_Aug18-2021.fa"
  # target.hyloidea = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Pipeline/Source_Files/Probe_Sets/Hyloidea_Aug18-2021.fa"
  # sample.groups = read.csv("sample_probeset_group.csv", header = T)
  # ranoidea.groups = sample.groups[sample.groups[,1] %in% "Ranoidea-V1",]
  # hyloidea.groups = sample.groups[sample.groups[,1] %in% "Hyloidea-V1",]
  # reduced.groups = sample.groups[sample.groups[,1] %in% "RedRan",]
  # mafft.path = "/usr/local/bin"
  # output.name = "reduced-genetic-distance"
  # target.type = "reference"
  #
  # threads = 4
  # memory = 4
  # overwrite = TRUE
  # quiet = FALSE
  #
  # alignment.directory = align.reduced
  # target.fasta = target.reduced


  ####################################################################
  ##### Required program path check
  ####################################################################
  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  ####################################################################
  ##### Input path check
  ####################################################################

  if (is.null(alignment.directory) == TRUE){ stop("Please provide input directory.") }

  if (target.type == "Reference" && is.null(target.fasta) == TRUE){ stop("Please provide reference fasta.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.name) == F){ dir.create(output.name) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.name))
      dir.create(output.name)
    }
  }#end else

  ####################################################################
  ##### Raw data collection for each alignment and sample
  ####################################################################

  #Gather alignments
  alignment.files = list.files(alignment.directory, full.names = T)
  alignment.files = alignment.files[grep(".phy$", alignment.files)]
  target.markers = Biostrings::readDNAStringSet(target.fasta)
  target.markers = target.markers[duplicated(names(target.markers)) != T]

  marker.names = gsub(".*\\/", "", alignment.files)
  marker.names = gsub(".phy$", "", marker.names)

  target.markers = target.markers[names(target.markers) %in% marker.names]
  marker.names = marker.names[marker.names %in% names(target.markers)]

  #Finds taxa
  taxa.temp = c()
  for (i in 1:(0.25*length(alignment.files))){
    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = alignment.files[i], format = "phylip"))
    taxa.temp = unique(append(taxa.temp, names(align)))
  }

  sample.names = unique(taxa.temp)

  #Sets up multiprocessing
  cl = parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  #all.data = c()
  all.data = foreach::foreach(i=1:length(marker.names),  .combine = rbind, .packages = c("PhyloCap", "foreach", "Biostrings","data.table", "ape", "stringr")) %dopar% {
  #Loops through each alignment
  #for (i in 1:length(marker.names)){

    #START HERE
    align.file = alignment.files[grep(paste0(marker.names[i], ".phy$"), alignment.files)]
    align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = align.file, format = "phylip"))

    if (target.type == "consensus"){

      #Make a consensus sequence
      con.seq = makeConsensus(alignment = align,
                              method = "majority",
                              type = "DNA")
      names(con.seq) = "Reference_Marker"
      new.align = append(align, con.seq)
    }#end target consensus

    if (target.type == "reference"){

      temp.target = target.markers[names(target.markers) == marker.names[i],]
      temp.target = temp.target[1]
      names(temp.target) = "Reference_Marker"

      #Aligns and then reverses back to correction orientation
      alignment = runMafft(sequence.data = align,
                           add.contigs = temp.target,
                           algorithm = "add",
                           adjust.direction = TRUE,
                           threads = 1,
                           cleanup.files = T,
                           quiet = quiet,
                           mafft.path = mafft.path)

      #Aligns and then reverses back to correction orientation
      reversed = names(alignment)[grep(pattern = "_R_", names(alignment))]
      if (length(reversed[grep(pattern = "Reference_Marker", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
      names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))

      #Removes the edge gaps
      ref.aligned = as.character(alignment['Reference_Marker'])
      not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
      ref.start = min(not.gaps)
      ref.finish = max(not.gaps)
      trim.align = Biostrings::subseq(alignment, ref.start, ref.finish)

      new.align = trimAlignmentColumns(alignment = trim.align,
                                       min.gap.percent = 50)

    }#end reference if

    #Finds the alignemnt pairwise distance from the target
    diff = pairwiseDistanceTarget(alignment = new.align,
                                  target = "Reference_Marker")

    sample.dist = diff[names(diff) != "Reference_Marker"]

    #Sets up data to collect
    header.all = c("sample" , "target_marker", "target_length", "genetic_distance")
    sample.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.all)))
    data.table::setnames(sample.data, header.all)
    sample.data[, sample:=as.character(sample.names)]
    sample.data[, target_marker:=as.character(target_marker)]
    sample.data[, target_marker:=as.character(marker.names[i])]
    sample.data[, target_length:=Biostrings::width(temp.target)]

    data.table::set(sample.data, i = match(names(sample.dist), sample.data$sample), j = match("genetic_distance", header.all), value = sample.dist )

    #Writes sample table
    #all.data = rbind(all.data, sample.data)
    print(data.frame(sample.data))

  }#end i loop

  parallel::stopCluster(cl)

  #Writes final table
  write.table(all.data, file = paste0(output.name, "/sample-genetic-distance_raw-data.txt"), sep = "\t", row.names = F)

  print(paste0("Saved raw data at ", output.name, "/sample-genetic-distance_raw-data.txt"))

  ####################################################################
  ##### Summarizes the raw data previously collected into a summary table
  ####################################################################

  #Get stats
  temp.min = unlist(lapply( split(all.data$genetic_distance, all.data$sample), min))
  temp.max = unlist(lapply( split(all.data$genetic_distance, all.data$sample), max))
  temp.mean = unlist(lapply( split(all.data$genetic_distance, all.data$sample), mean))
  temp.median = unlist(lapply( split(all.data$genetic_distance, all.data$sample), median))
  temp.sd = unlist(lapply( split(all.data$genetic_distance, all.data$sample), sd))

  #Sets up data to collect
  header.summ = c("sample", "mean_genetic_distance", "median_genetic_distance",
                 "min_genetic_distance", "max_genetic_distance", "sd_genetic_distance")
  summary.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.summ)))
  data.table::setnames(summary.data, header.summ)
  summary.data[, sample:=as.character(sample)]
  summary.data[, sample:=as.character(sample.names)]

  data.table::set(summary.data, i = match(names(temp.min), summary.data$sample), j = match("min_genetic_distance", header.summ), value = temp.min )
  data.table::set(summary.data, i = match(names(temp.max), summary.data$sample), j = match("max_genetic_distance", header.summ), value = temp.max )
  data.table::set(summary.data, i = match(names(temp.mean), summary.data$sample), j = match("mean_genetic_distance", header.summ), value = temp.mean )
  data.table::set(summary.data, i = match(names(temp.median), summary.data$sample), j = match("median_genetic_distance", header.summ), value = temp.median )
  data.table::set(summary.data, i = match(names(temp.sd), summary.data$sample), j = match("sd_genetic_distance", header.summ), value = temp.sd )

  write.table(summary.data, file = paste0(output.name, "/sample-genetic-distance_summary.txt"), sep = "\t", row.names = F)

  print(paste0("Saved raw data at ", output.name, "/sample-genetic-distance_summary.txt"))

}#end function

### END SCRIPT
