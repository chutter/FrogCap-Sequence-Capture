import nimib, nimibook

nbInit

nbUseNimibook 

nbText: md"""
# Data Processing
This tutorial will guide you through the processing of your raw data into assembled contigs for each of your samples. 

## Create simple configuration file

The configuration file is quite simple for the FrogCap pipeline to function properly, and can even be excluded. Fastp trims the raw reads of adapters without needing a file of adapter sequences or barcodes, as it can detect adapter contamination algorithmically. 

This is the "File_rename.csv" file, which has only two columns: File and Sample.

The "File" column: the unique string that is part of the file name for the two read pairs, while excluding read information. Example:

> ``Spinomantis_elegans_CRH111_AX1212_L001_R1.fastq.gz``

> ``Spinomantis_elegans_CRH111_AX1212_L001_R2.fastq.gz``

Are the two sets of reads for a given sample. Your "File" column value would then be:

> ``Spinomantis_elegans_CRH111``(or Spinomantis_elegans_CRH111_AX1212; both will work)


The "Sample" column: What you would like your sample name to be. This will be used up to alignments and trees. Ensure that your samples all have unique names and are not contained within each other (e.g. Genus_species_0, Genus_species_01). Also exclude special characters and replace spaces with underscores. Hyphens are also ok. 
In this example, the "Sample" Column would be: 
>
>Spinomantis_elegans_CRH111
>

All put together:


File | Sample
------|:--------:|
Spinomantis_elegans_CRH111 | Spinomantis_elegans_CRH111
Aglyptodactylus_securifer_CRH1644 | Aglyptodactylus_securifer_CRH1644
Boophis_burgeri_CRH0481 | Boophis_burgeri_CRH0481
Mantidactylus_femoralis_CRH2340 | Mantidactylus_femoralis_CRH2340

Note that the columns will not be the same if the sample name is not contained in your file name. This may look like: 

> ``CRH111_AX1212_L001_R1.fastq.gz``

> ``CRH111_AX1212_L001_R1.fastq.gz``



File | Sample
------|:--------:|
CRH111_AX1212 | Spinomantis_elegans_CRH111
CRH1644_AX1212 | Aglyptodactylus_securifer_CRH1644
CRH0481_AX1212 | Boophis_burgeri_CRH0481
CRH2340_AX1212 | Mantidactylus_femoralis_CRH2340


## Raw data input 

Finally, before running the FrogCap pipeline scripts, the input files must be organized in a specific way. First, create a directory and name it after the project name. Second, put the newly created "File_rename.csv", the downloaded and unzipped folder of "FrogCap_Files", and finally the demultiplexed raw reads in a folder named "raw_data". An example of the file structure is shown below.


```
     /Project_Name
      ├── /raw_data
      │   ├── Spinomantis_elegans_CRH111_R1.fastq.gz
      │   ├── Spinomantis_elegans_CRH111_R2.fastq.gz
      │   ├── Boophis_burgeri_CRH0481_R1.fastq.gz
      │   ├── Boophis_burgeri_CRH0481_R2.fastq.gz
      │   ├── Aglyptodactylus_securifer_CRH1644_R1.fastq.gz
      │   ├── Aglyptodactylus_securifer_CRH1644_R2.fastq.gz
      │   ├── Mantidactylus_femoralis_CRH2340_R1.fastq.gz
      │   └── Mantidactylus_femoralis_CRH2340_R2.fastq.gz
      ├── File_rename.csv
      └── /FrogCap_Files

*/ denotes directory

```


## Preprocess the Reads [R script 01_Pre_Process_Reads.R]

The first script that will be used is "01_Pre_Process_Reads.R". Open this script in any text or R editor to edit the run parameters. 

```

##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#Cluster
work.dir<-"/home/Project_Name"
raw.dir<-"/home/Project_Name/raw_reads"
decontam.path<-"/home/Project_Name/FrogCap_Files/Contamination_Genomes" #folder of contaminants I gave

#Various settings of choice
#Sample file
# "File" column for you file name
# "Sample" column should correspond to the file name for each samples set of reads e.g. "CRH01"
sample.file<-"File_rename.csv" #two column file with "File" and "Sample"
out.dir<-"Processed_Samples" #I would keep this name. 
threads<-10 #number of threads, 8-10 is a recommended amount
mem<-"80g" #need the "g' character here

```

This section is the only section that needs editing. Namely, the directories should be the full path to make things go smoother. Also note the threads and memory that can be adjusted. 


When finished editing, run the script using the "Rscript" terminal command which installs alongside R. This can be done manually in the terminal window or added as the command in a cluster file. 

```

> Rscript 01_Pre_Process_Reads.R

```

This script takes about a 30-45 minutes a sample using the threads and memory from above. When completed, it should create several types of new files and directories, as shown here:

```
     /Project_Name
      ├── /raw_data
      ├── File_rename.csv
      ├── /FrogCap_Files
      └── /Processed_Samples
          ├── /Spinomantis_elegans_CRH111
          ├── /Boophis_burgeri_CRH0481
          ├── /Aglyptodactylus_securifer_CRH1644
          └── /Mantidactylus_femoralis_CRH2340

*/ denotes directory

```

Inside the Processed_Samples directory is a directory for every sample. In each of these sample directory are the processed reads and each sample should have the following files: 


```

     /Sample_Name_ID
      ├── /reports
      │   ├── read_merging_stats.csv
      │   ├── cleaned-reads-snp_stats.csv
      │   ├── Adapter_trimming_fastp.html
      │   └── Adapter_trimming_fastp.json
      ├── /assembly-reads
      │   ├── Sample_name_ID_READ1.fastq.gz
      │   ├── Sample_name_ID_READ2.fastq.gz
      │   └── Sample_name_ID_singletons.fastq.gz
      ├── /cleaned-reads-snp
      │   ├── Sample_name_ID_READ1.fastq.gz
      │   └── Sample_name_ID_READ2.fastq.gz
      └── /raw-reads-symlink
          ├── Sample_name_ID_READ1.fastq.gz
          └── Sample_name_ID_READ2.fastq.gz

*/ denotes directory

```

The directories contains reads filtered differently for different types of analyses:

* assembly-reads: these reads are for SPades to assemble. They have been cleaned, filtered, and merged. 

* cleaned-reads-snp: these reads are for downstream SNP calling analyses. They have only been cleaned. 

* raw-reads-symlink: these are symlinks (i.e. shortcuts) to your raw reads and don't use up much space. 

* QC-files: Summary files from processing steps, and other summary statistics will be placed here. 


## Assemble the processed reads [R script 02_Assemble_Spades.R]

After cleaning and filtering from the previous step, the second script will assemble the cleaned reads. The name of the script that will be used is "02_Assemble_Spades.R". Open this script in any text or R editor to edit the run parameters. 


```r
##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#Set up parameters
threads<-"8"
mem<-"80"

#Directory where the "Processed_Samples" Folder is located
work.dir<-"/home/Project_Name"

######################## Probably shouldn't change these 
#Choose which directory to use for assembly:
#aseembly-reads = recommended assembly reads directory (default)
#cleaned-reads-snp = if the assembly_reads are too aggressively filtered. 
read.dir<-"assembly-reads"
out.dir<-"Assembled_Contigs"
in.dir<-"Processed_Samples"
clean<-"True" #Deletes intermediate files, recommended
```

In the "Probably shouldn't change these" section, changing the assembly directory from the "assembly-reads" to the "cleaned-reads-snp" can help in situations where the "assembly-reads" give poor results. This can happen when the raw data is low depth of coverage and quality. The other folders are existing directory names if they have been changed from the defaults. "Clean" will remove intermediate assembly files, only change if you need a specific output file. 

When finished editing, run the script using the "Rscript" terminal command which installs alongside R. This can be done manually in the terminal window or added as the command in a cluster file. 

```bash
> Rscript 02_Assemble_Spades.R
```

This script can take up to 10 hours a sample, but averages around 2-3 (using the threads and memory from above). Its also possible to run multiple instances of the script on a subset of the samples to speed things along. When completed, it should create three new files and one new directory, as shown here:

```

     /Sample_Name_ID
      ├── /reports
      │   ├── assembly-reads_stats.csv *** NEW
      │   ├── read_merging_stats.csv
      │   ├── cleaned-reads-snp_stats.csv
      │   ├── Adapter_trimming_fastp.html
      │   └── Adapter_trimming_fastp.json
      ├── /assembled-contigs *** NEW
      │   ├── Sample_name_ID_contigs.fa *** NEW
      │   └── Sample_name_ID_dipcontigs.fa *** NEW
      ├── /assembly-reads
      │   ├── Sample_name_ID_READ1.fastq.gz
      │   ├── Sample_name_ID_READ2.fastq.gz
      │   └── Sample_name_ID_singletons.fastq.gz
      ├── /cleaned-reads-snp
      │   ├── Sample_name_ID_READ1.fastq.gz
      │   └── Sample_name_ID_READ2.fastq.gz
      └── /raw-reads-symlink
          ├── Sample_name_ID_READ1.fastq.gz
          └── Sample_name_ID_READ2.fastq.gz

*/ denotes directory

```

The Assembly run generated a statistics file in the "reports" folder and a new folder "assembled-contigs". This folder has two new files: the "_contigs.fa" file are the Spades assembled contigs while the "_dipcontigs.fa" file are the regular contigs where contigs from the same haplotype are merged together to create a "haplocontig". 

These two scripts can be run on any number of samples in any combination. However, when done with Tutorial 1, you will need to select your definitive sample set for alignment for Tutorial 2. 

"""

nbSave