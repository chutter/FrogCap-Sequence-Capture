#' @title summary.fastqStats
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

summary.fastqStats = function(read.directory = NULL,
                      sub.directory = NULL,
                      output.name = "fastq-stats",
                      read.length = 150,
                      threads = 1,
                      mem = 1,
                      overwrite = FALSE) {

  # #Debug
  #  setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
  # read.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Processed_Samples"
  #  sub.directory = "assembly-reads"
  #  output.name = "fastq-stats"
  #  read.length = 150
  #  overwrite = TRUE

  #Quick checks
  if (is.null(read.directory) == TRUE){ stop("Please provide input reads.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.name) == F){ dir.create(output.name) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.name))
      dir.create(output.name)
    }
  }#end else

  #Read in sample data **** sample is run twice?!
  reads = list.files(read.directory, recursive = T, full.names = T)
  if (is.null(sub.directory) != TRUE) {
    reads = reads[grep(paste0(sub.directory, "/"), reads)]
    sample.names = gsub(paste0("/", sub.directory, "/.*"), "", reads)
    sample.names = unique(gsub(paste0(read.directory, "/"), "", sample.names))
  } else {
    sample.names = list.files(read.directory, recursive = F, full.names = F)
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }


  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Read1_Count = as.numeric(),
                             Read2_Count = as.numeric(),
                             Read3_Count = as.numeric(),
                             Total_Reads = as.numeric(),
                             Read_Pairs = as.numeric(),
                             Read_Length = as.numeric(),
                             MegaBasePairs = as.numeric(),
                             Reads_Per_Million = as.numeric())

  #Runs through each sample
  for (i in 1:length(sample.names)) {
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

      #Gathers stats on initial data
      read1.count = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
      read2.count = as.numeric(system(paste0("zcat < ", lane.reads[2], " | echo $((`wc -l`/4))"), intern = T))
      if (length(lane.reads) == 3){
        read3.count = as.numeric(system(paste0("zcat < ", lane.reads[3], " | echo $((`wc -l`/4))"), intern = T))
      } else { read3.count = 0 }

      scale.factor = (read1.count + read2.count + read3.count) / 1000000
      if (read1.count == read2.count){ read.pairs = read1.count } else { read.pairs = "PairsUnequal"}

      temp.remove = data.frame(Sample = sample.names[i],
                               Read1_Count = read1.count,
                               Read2_Count = read2.count,
                               Read3_Count = read3.count,
                               Total_Reads = read1.count + read2.count,
                               Read_Pairs = read.pairs,
                               Read_Length = read.length,
                               MegaBasePairs = (read.length * (read1.count + read2.count + read3.count) )/1000000,
                               Reads_Per_Million = scale.factor)

      summary.data = rbind(summary.data, temp.remove)


    }#end sample j loop

    print(paste0(sample.names[i], " Completed fastq counting!"))

  }#end sample i loop

  write.csv(summary.data, file = paste0(output.name, ".csv"), row.names = FALSE)
  return(summary.data)
}

