#' @title readDepth
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

readDepth = function(read.directory = NULL,
                     sub.directory = NULL,
                     output.name = "read-depth",
                     output.binned = FALSE,
                     threads = 1,
                     memory = 1,
                     overwrite = FALSE,
                     quiet = TRUE,
                     samtools.path = NULL,
                     bwa.path = NULL,
                     picard.path = NULL) {

  #Debug
  # setwd("/Volumes/Armored/FrogCap_Anura_Seqcap/Analyses")
  # read.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Processed_Samples"
  # sub.directory = "cleaned-reads-snp"
  # output.name = "read-depth"
  # overwrite = FALSE
  # samtools.path = "/Users/chutter/miniconda3/bin"
  # bwa.path = "/usr/local/bin"
  # picard.path = "/Users/chutter/miniconda3/bin"
  # bbmap.path = "/usr/local/bin"
  # threads = 4
  # memory = 4
  # output.binned = FALSE
  # quiet = TRUE

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
      system(paste0("rm ", output.name))
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

  #Sets up data to collect
  header.data = c("Sample","numberLoci", "meanDepth", "medianDepth", "sdDepth", "minDepth", "maxDepth")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Sample:=as.character(Sample)]

  #Sets up bin data
  bin.collect = data.table::data.table(matrix(as.numeric(0), nrow = 100, ncol = length(sample.names)))
  data.table::setnames(bin.collect, gsub(".*\\/", "", sample.names))

  #Loops through each locus and does operations on them
  for (i in 1:length(sample.names)){

    #Reads in the contig file to check
    dir.create(paste0(output.name, "/", sample.names[i]))
    sample.contigs = paste0(read.directory, "/", sample.names[i], "/assembled-contigs/", sample.names[i], "_orthologs.fa")

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
      system(paste0("cp ", sample.contigs, " ", output.name, "/", sample.names[i], "/sample-reference.fa"))
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
                    " BuildBamIndex INPUT=", output.name, "/", sample.names[i], "/dedup_paired.bam",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #Gets the genome coverage using bedtools
      # NOTE CAN USE -a OPTION TO GET ALL POSITIONS
      system(paste0(samtools.path, "samtools depth -aa -m -s ",
                    output.name, "/", sample.names[i], "/dedup_paired.bam > ",
                    output.name, "/", sample.names[i], "/samtools_perbase_depth.txt"))
      system(paste0(samtools.path, "samtools idxstats --threads ", threads, " ",
                    output.name, "/", sample.names[i], "/dedup_paired.bam > ",
                    output.name, "/", sample.names[i], "/samtools_idxstats.txt"))

      system(paste0("rm ", output.name, "/", sample.names[i], "/paired.bam"))

      #Loads in depth information
      d.headers = c("locus", "base", "mapped_reads")
      depth.stats = data.table::fread(file = paste0(output.name, "/", sample.names[i], "/samtools_perbase_depth.txt"))
      data.table::setnames(depth.stats, d.headers)

      i.headers = c("locus", "length", "mapped_reads", "unmapped_reads")
      idx.stats = data.table::fread(file = paste0(output.name, "/", sample.names[i], "/samtools_idxstats.txt"))
      data.table::setnames(idx.stats, i.headers)
      idx.stats = idx.stats[idx.stats$locus != "*",]

      #Goes through each locus and calculates stats
      locus.names = unique(depth.stats$locus)
      total.mapped.reads = sum(idx.stats$mapped_reads)

      #Starter stats
      #Sets up data to collect
      header.spp = c("sample", "locus", "locus_length", "mapped_reads", "locus_rpkm",
                     "mean_depth", "median_depth", "sd_depth", "min_depth", "max_depth")
      spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(locus.names), ncol = length(header.spp)))
      data.table::setnames(spp.data, header.spp)
      spp.data[, sample:=as.character(sample)]
      spp.data[, locus:=as.character(locus)]

      #summary stats
      locus.rpkm = idx.stats$mapped_reads / (idx.stats$length/1000 * total.mapped.reads/1000000)
      temp.mean = unlist(lapply( split(depth.stats$mapped_reads, depth.stats$locus), mean))
      temp.median = unlist(lapply( split(depth.stats$mapped_reads, depth.stats$locus), median))
      temp.sd = unlist(lapply( split(depth.stats$mapped_reads, depth.stats$locus), sd))
      temp.min = unlist(lapply( split(depth.stats$mapped_reads, depth.stats$locus), min))
      temp.max = unlist(lapply( split(depth.stats$mapped_reads, depth.stats$locus), max))

      #Main stasts
      data.table::set(spp.data, j = match("sample", header.spp), value = sample.names[i] )
      data.table::set(spp.data, j = match("locus", header.spp), value = locus.names )
      data.table::set(spp.data, i = match(locus.names, idx.stats$locus), j = match("locus_length", header.spp), value = idx.stats$length )
      data.table::set(spp.data, i = match(locus.names, idx.stats$locus), j = match("mapped_reads", header.spp), value = idx.stats$mapped_reads )
      data.table::set(spp.data, i = match(locus.names, idx.stats$locus), j = match("locus_rpkm", header.spp), value = locus.rpkm )

      #Summary stats
      data.table::set(spp.data, i = match(locus.names, names(temp.mean)), j = match("mean_depth", header.spp), value =  temp.mean)
      data.table::set(spp.data, i = match(locus.names, names(temp.median)), j = match("median_depth", header.spp), value =  temp.median)
      data.table::set(spp.data, i = match(locus.names, names(temp.sd)), j = match("sd_depth", header.spp), value =  temp.sd)
      data.table::set(spp.data, i = match(locus.names, names(temp.min)), j = match("min_depth", header.spp), value =  temp.min)
      data.table::set(spp.data, i = match(locus.names, names(temp.max)), j = match("max_depth", header.spp), value =  temp.max)

      #Save the two datasets inside the depth_analysis folder
      write.table(spp.data, file = paste0(output.name, "/", sample.names[i], "/species_summary_data.txt"), sep = "\t", row.names = F)

      #Confirm completed
      print(paste0(sample.names[i], " Completed!"))

    }#end j loop

  }#end i loop

 # write.table(collect.data, file = "Sample_Coverage_Summary.txt", sep = "\t", row.names = F)

}#end function
