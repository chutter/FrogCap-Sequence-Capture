#' @title summarizeAlignments
#'
#' @description Function for gather summary statistics on your alignment. Can be used for filtering or summarizing data.
#'
#' @param alignment.path path to a folder of sequence alignments in phylip format.
#'
#' @param file.export give a save name if you wnat to save the summary to file.
#'
#' @param overwrite if TRUE overwrites file if it exists; FALSE the dataset is skipped
#'
#' @param dataset.name A name for your dataset. i.e. exons, introns, UCEs
#'
#' @param alignment.type select the format of the alignment. Phylip is avaialble for now, will be expanded in the future.
#'
#' @return returns a data.table with the raw summary statistics calculated for each alignment in the alignment.path. A csv file can optionally be saved by giving a file name to file.export
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


summarizeAlignments = function(alignment.path = NULL,
                               file.export = NULL,
                               overwrite = FALSE,
                               dataset.name = NULL,
                               alignment.format = c("phylip", "nexus")) {
  #Parameter checks
  if(is.null(alignment.path) == TRUE){ stop("Error: no alignment path provided.") }
  if(is.null(dataset.name) == TRUE){ stop("Error: a dataset name is needed.") }

  #Check if files exist or not
  if (dir.exists(alignment.path) == F){
    return(paste0("Directory of alignments could not be found. Exiting."))
  }#end file check

  #Overwrite checker
  if (overwrite == TRUE){
    if (file.exists(paste0(file.export, ".csv")) == T){
      #Checks for output directory and creates it if not found
      system(paste0("rm ", file.export, ".csv"))
    }#end file exists
  } else {
    if (file.exists(paste0(file.export, ".csv")) == T){
      print(paste0("File exists for ", file.export, " and overwrite = FALSE. Exiting."))
      save.data = read.csv(paste0(file.export, ".csv"))
      return(save.data)
    }#end file check
  }#end else

  #Gets list of alignments from path
  align.names = list.files(alignment.path)

  #Collects the super cool data
  header.data = c("dataset", "file", "number_samples", "proportion_samples", "alignment_length",
                  "count_pis", "proportion_pis", "count_missing_bp", "proportion_missing_bp")
  #Sets up data collection data.frame
  collect.data = data.table::data.table(matrix(as.numeric(0),
                                   nrow = length(align.names),
                                   ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, file:=as.character(file)]
  collect.data[, dataset:=as.character(dataset)]

  #Loops through each alignment to gather statistics
  for (x in 1:length(align.names)){
    #Reads in alignment

    if (alignment.format == "phylip"){
      align = ape::read.dna(paste0(alignment.path, "/", align.names[x]),
                       format = "sequential")
    }
    if (alignment.format == "nexus"){
      align = ape::read.nexus.data(paste0(alignment.path, "/", align.names[x]))
      align = ape::as.DNAbin(matrix(unlist(align), ncol = length(align[[1]]), byrow = TRUE))
    }

    #Collect data
    set(collect.data, i = as.integer(x), j = match("dataset", header.data), value = dataset.name )
    set(collect.data, i = as.integer(x), j = match("file", header.data), value = align.names[x] )
    #Sample data
    set(collect.data, i = as.integer(x), j = match("number_samples", header.data), value = nrow(align) )
    #Length data
    set(collect.data, i = as.integer(x), j = match("alignment_length", header.data), value = ncol(align) )

    count.pis = informativeSites(align, count = T, ambiguities = T)
    prop.pis = round(count.pis/ncol(align),3)
    set(collect.data, i = as.integer(x), j = match("count_pis", header.data), value = count.pis)
    set(collect.data, i = as.integer(x), j = match("proportion_pis", header.data), value = prop.pis)

    #Removes samples that too short individually
    len.temp = as.character(as.list(align))
    len.loci = lapply(len.temp, function (x) x[x != "-"])
    len.loci = lapply(len.loci, function (x) x[x != "n"])
    len.loci = lapply(len.loci, function (x) x[x != "?"])
    spp.len = unlist(lapply(len.loci, function (x) length(x)))
    miss.total = (max(spp.len) - spp.len)
    miss.prop = round(sum(miss.total)/(max(spp.len)*nrow(align)), 3)

    #Get missing bp data
    set(collect.data, i = as.integer(x), j = match("count_missing_bp", header.data), value =  sum(miss.total) )
    set(collect.data, i = as.integer(x), j = match("proportion_missing_bp", header.data), value = miss.prop )

  } #x loop

  save.data = collect.data[collect.data$file != 0,]
  save.data[, proportion_samples:=round(number_samples/max(number_samples), 3)]

  if (is.null(file.export) != TRUE){
    write.csv(save.data, file = paste0(file.export, ".csv"), row.names = F)
  }#end if

  return(save.data)

}#End function summarizeAlignments

