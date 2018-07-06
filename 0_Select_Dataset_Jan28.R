###### This script creates a new project folder if you have lots of data in many different folders within a single directory

#HOME 
proj.name<-"Microhylidae_Dataset"
#Directory locations
raw.dir<-"/Volumes/Armored/Raw_Data"
setwd(raw.dir)

read.files<-list.files(".", recursive = T)
read.files<-read.files[grep(pattern = ".fastq.gz$", x=read.files)]

#Various settings of choice
sample.file<-"/Volumes/Armored/Raw_Data/Datasheets/Master_Datasets_Feb9.csv"

#Creates output directory
dir.create(proj.name)

#Reads in the dataset spreadsheet to choose which to save
sample.data<-read.csv(sample.file, stringsAsFactors = F)
copy.data<-sample.data[colnames(sample.data)%in% proj.name] 
sub.data<-sample.data[copy.data[,1] == 1,]

#Copy and rename files to output directory
for (i in 1:nrow(sub.data)){
  
  sample.reads<-read.files[grep(pattern = sub.data$File[i], x = read.files)]
  
  if (length(sample.reads) == 1){ 
    print(paste("only 1 read for ", sub.data$Sample[i], sep = ""))
    next
    }
  
  folder.pick<-unique(gsub(pattern= "/.*", replacement = "", x = sample.reads))[1]
  
  copy.reads<-sample.reads[grep(pattern = folder.pick, x= sample.reads)]
  
  system(paste("cp ", copy.reads[1], " ", raw.dir, "/", proj.name, "/", sub.data$Sample[i], "_R1.fastq.gz", sep = ""))
  system(paste("cp ", copy.reads[2], " ", raw.dir, "/", proj.name, "/", sub.data$Sample[i], "_R2.fastq.gz", sep = ""))

}
