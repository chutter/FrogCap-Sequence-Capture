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

binnedDepthSummary = function(depth.directory = NULL,
                              output.directory = "binned-depth-summary",
                              sample.groups = NULL,
                              threads = 1,
                              memory = 1,
                              overwrite = FALSE) {

  #Debug
  setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
  depth.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Analyses/read-depth"
  sample.groups = "sample_probeset_group.csv"
  output.directory = "Analyses/binned-depth-summary"
  overwrite = TRUE
  max.value = 40000000

  library(ggplot2)
  parameters = theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(vjust = 1, hjust = 0.5, size = 28)) +
    theme(axis.title.x = element_text(size = 20, vjust=-0.5),
          axis.title.y = element_text(size = 20, vjust=2.2)) +
    theme(text = element_text(size=20),
          axis.text.x = element_text(vjust=1))


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

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  group.data = read.csv(sample.groups, header = T)
  group.data = group.data[group.data[,2] %in% sample.names,]
  group.names = unique(group.data[,1])

  ##################################### Sample Summary stats
  ###############################################################################################
  #Sets up data to collect
  header.spp = c("group", "sample", "mean_rpkm", "median_rpkm", "sd_rpkm", "min_rpkm", "max_rpkm")
  spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.spp)))
  data.table::setnames(spp.data, header.spp)
  spp.data[, sample:=as.character(sample)]
  spp.data[, group:=as.character(group)]

  for (j in 1:length(group.names)){

    #Runs through each sample
    #all.data = data.table::data.table()

    group.samples = group.data[group.data[,1] %in% group.names[j],][,2]
    bin.vals = vector("list", 100)
    for (i in 1:length(group.samples)) {
      #################################################
      ### Part A: prepare for loading and checks
      #################################################
      if (is.null(sub.directory) == TRUE){
        sample.data = as.matrix( data.table::fread(paste0(depth.directory, "/", group.samples[i], "/rpkm-binned-data.txt")) )
      } else {
        sample.data = as.matrix( data.table::fread(paste0(depth.directory, "/",group.samples[i], "/", sub.directory, "/rpkm-binned-data.txt")) )
      }

      sample.median = apply(sample.data, 1, median, na.rm=T)
      bin.vals = mapply(c, bin.vals, sample.median, SIMPLIFY=FALSE)

    } #end i loop

    #Gets barplot data.frame ready
    bin.data = data.frame(bin = seq(1:100), median = unlist(lapply(bin.vals, median)),
                         sd = unlist(lapply(bin.vals, sd)) )

    #ggplot rocks
    ggplot(bin.data, aes(x=bin, y=median)) + parameters +
      geom_bar(position=position_dodge(), stat="identity", size=.3, color = "#00FF7F", fill = "#00FF7F") +
      geom_errorbar(aes(ymin=median-sd, ymax=median+sd), size=.1, width=.1, position=position_dodge(.9), color = "black") +
      ggtitle(group.names[j]) + xlab("Contig (1 percent bins)") + ylab("Median Depth (RPKM)") + ylim(0, max.value)

    ggplot2::ggsave(paste0(output.directory, "/", group.names[j], "_bin-plot_rpkm-median.pdf"), width = 8, height = 6)

  }#end j loop


}#end function





