#' @title summarizeDepth
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

summarizeDepth = function(depth.directory = NULL,
                          sub.directory = NULL,
                          output.directory = "depth-summary",
                          sample.groups = NULL,
                          remove.samples = NULL,
                          threads = 1,
                          mem = 1,
                          overwrite = FALSE) {

  #Debug
  setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
  depth.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Analyses/read-depth"
  sub.directory = NULL
  sample.groups = "sample_probeset_group.csv"
  output.directory = "Analyses/depth-summary"
  remove.samples = c("Leptobrachium_hainensis_KUFS_194", "Bolitaglossa_palmata_LAC_1154")
  overwrite = TRUE

  #Quick checks
  if (is.null(depth.directory) == TRUE){ stop("Please provide input directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Read in sample data **** sample is run twice?!
  stat.files = list.files(depth.directory, recursive = T, full.names = T)
  if (is.null(sub.directory) == FALSE) {
    stat.files = stat.files[grep(paste0(sub.directory, "/"), stat.files)]
    sample.names = gsub(paste0("/", sub.directory, "/.*"), "", stat.files)
    sample.names = unique(gsub(paste0(depth.directory, "/"), "", sample.names))
  } else {
    sample.names = list.files(depth.directory, recursive = F, full.names = F)
  }

  sample.names = sample.names[!sample.names %in% remove.samples]

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  group.data = read.csv(sample.groups, header = T)

  ##################################### Sample Summary stats
  ###############################################################################################
  #Sets up data to collect
  header.spp = c("group", "sample", "mean_rpkm", "median_rpkm", "sd_rpkm", "min_rpkm", "max_rpkm")
  spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.spp)))
  data.table::setnames(spp.data, header.spp)
  spp.data[, sample:=as.character(sample)]
  spp.data[, group:=as.character(group)]

  #Runs through each sample
  all.data = data.table::data.table()
  for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    if (is.null(sub.directory) == TRUE){
      sample.data = data.table::fread(paste0(depth.directory, "/", sample.names[i], "/species_summary_data.txt"))
    } else {
      sample.data = data.table::fread(paste0(depth.directory, "/",sample.names[i], "/", sub.directory, "/species_summary_data.txt"))
    }

    #Main stasts
    temp.group = group.data[group.data$Sample == sample.names[i],]
    data.table::set(spp.data, i = as.integer(i), j = match("group", header.spp), value = temp.group[1,1] )
    data.table::set(spp.data, i = as.integer(i), j = match("sample", header.spp), value = sample.names[i] )

    #Summary stats
    data.table::set(spp.data, i = as.integer(i), j = match("mean_rpkm", header.spp), value =  mean(sample.data$locus_rpkm) )
    data.table::set(spp.data, i = as.integer(i), j = match("median_rpkm", header.spp), value =  median(sample.data$locus_rpkm) )
    data.table::set(spp.data, i = as.integer(i), j = match("sd_rpkm", header.spp), value =  sd(sample.data$locus_rpkm) )
    data.table::set(spp.data, i = as.integer(i), j = match("min_rpkm", header.spp), value =  min(sample.data$locus_rpkm) )
    data.table::set(spp.data, i = as.integer(i), j = match("max_rpkm", header.spp), value =  max(sample.data$locus_rpkm) )

    sample.data[, group:=as.character(temp.group[1,1])]

    all.data = rbind(all.data, sample.data)

  }#end i loop

  write.csv(spp.data, file = paste0(output.directory, "/sample_rpkm-summary.csv"), row.names = FALSE, quote = FALSE)

  ################################################ Comparison plots

  ggplot2::ggplot(data = spp.data, ggplot2::aes(x =group, y = mean_rpkm)) +
    ggplot2::geom_dotplot(ggplot2::aes(fill = group), binaxis='y', stackdir='centerwhole',
                          stackratio=0.9, dotsize=1.1) + ggplot2::theme_bw()

  ggplot2::ggsave(paste0(output.directory, "/DotPlot_Sample_rpkm-mean.pdf"), width = 9, height = 6)

  ggplot2::ggplot(data = spp.data, ggplot2::aes(x =group, y = median_rpkm)) +
    ggplot2::geom_dotplot(ggplot2::aes(fill = group), binaxis='y', stackdir='centerwhole',
                          stackratio=0.9, dotsize=1.2) + ggplot2::theme_bw()

  ggplot2::ggsave(paste0(output.directory, "/DotPlot_Sample_rpkm-median.pdf"), width = 9, height = 6)

  ggplot2::ggplot(data = spp.data, ggplot2::aes(x =group, y = median_rpkm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), width = 0.8) + ggplot2::theme_bw()


  ################################## Marker Summary stats
  ###############################################################################################

  group.names = unique(all.data$group )

  marker.data = data.table::data.table()
  for (k in 1:length(group.names)){

    #Gets group subset
    group.stats = all.data[all.data$group %in% group.names[k],]
    group.stats$locus = gsub("_\\|_.*", "", group.stats$locus)
    locus.names = unique(group.stats$locus)

    header.spp = c("group", "locus", "mean_length", "mean_rpkm", "median_rpkm", "sd_rpkm", "min_rpkm", "max_rpkm")
    spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(locus.names), ncol = length(header.spp)))
    data.table::setnames(spp.data, header.spp)
    spp.data[, group:=as.character(group)]
    spp.data[, locus:=as.character(locus)]

    #summary stats
    temp.mean = unlist(lapply( split(group.stats$locus_rpkm, group.stats$locus), mean))
    temp.median = unlist(lapply( split(group.stats$locus_rpkm, group.stats$locus), median))
    temp.sd = unlist(lapply( split(group.stats$locus_rpkm, group.stats$locus), sd))
    temp.min = unlist(lapply( split(group.stats$locus_rpkm, group.stats$locus), min))
    temp.max = unlist(lapply( split(group.stats$locus_rpkm, group.stats$locus), max))

    #Main stasts
    data.table::set(spp.data, j = match("group", header.spp), value = group.names[k] )
    data.table::set(spp.data, j = match("locus", header.spp), value = locus.names )

    #Summary stats
    data.table::set(spp.data, i = match(locus.names, names(temp.mean)), j = match("mean_rpkm", header.spp), value =  temp.mean)
    data.table::set(spp.data, i = match(locus.names, names(temp.median)), j = match("median_rpkm", header.spp), value =  temp.median)
    data.table::set(spp.data, i = match(locus.names, names(temp.sd)), j = match("sd_rpkm", header.spp), value =  temp.sd)
    data.table::set(spp.data, i = match(locus.names, names(temp.min)), j = match("min_rpkm", header.spp), value =  temp.min)
    data.table::set(spp.data, i = match(locus.names, names(temp.max)), j = match("max_rpkm", header.spp), value =  temp.max)

    marker.data = rbind(marker.data, spp.data)


  }#end group names j

  write.csv(marker.data, file = paste0(output.directory, "/marker_rpkm-summary.csv"), row.names = FALSE, quote = FALSE)

  ################################################ Comparison plots

  marker.data = marker.data[marker.data$median_rpkm < 500,]
  ggplot2::ggplot(data = marker.data, ggplot2::aes(x =group, y = median_rpkm)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), width = 0.8) + ggplot2::theme_bw()

  ggplot2::ggsave(paste0(output.directory, "/BoxPlot_Marker_rpkm-median.pdf"), width = 9, height = 6)


} #end function

