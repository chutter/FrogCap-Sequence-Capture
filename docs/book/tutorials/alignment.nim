import nimib, nimibook

nbInit

nbUseNimibook 

nbText: md"""
# Data Alignment

*** STILL A WORK IN PROGRESS

This tutorial will guide you through matching assembled contigs to the target markers and next aligning each marker. 

## Target Marker Matching [R script 03_Probe_Matching.R]

After assembly, confirm that each sample folder has a new "assembled-contigs" folder present. In this folder are two files: the "_contigs.fa" file are the Spades assembled contigs while the "_dipcontigs.fa" file are the regular contigs where contigs from the same haplotype are merged together to create a "haplocontig" from the Dipspades function of Spades. The dipcontigs are typically used in this pipeline, as it prevents data from being discarded for having the appearance of a paralog. However, the pre-merged contigs are available also. 

```
     /Spinomantis_elegans_CRH111
      ├── /reports
      ├── /assembled-contigs *** NEW
      │   ├── Spinomantis_elegans_CRH111_contigs.fa *** NEW
      │   └── Spinomantis_elegans_CRH111_dipcontigs.fa *** NEW
      ├── /assembly-reads
      ├── /cleaned-reads-snp
      └── /raw-reads-symlink

*/ denotes directory
```

Once you confirm that the assembly finished properly for each sample, check the "Project_Name" directory to check that the "Assembled_Contigs" folder is present and contains all of the samples that were assembled in the previous step. This is also the point where additional contigs from other assembly runs can be added in, as the sample contigs placed in the "Assembled_Contigs" folder will be the samples used for sequence alignment. 

```
     /Project_Name
      ├── /Assembled_Contigs *** NEW
      │   ├── Spinomantis_elegans_CRH111.fa *** NEW
      │   ├── Boophis_burgeri_CRH0481.fa *** NEW
      │   ├── Aglyptodactylus_securifer_CRH1644.fa *** NEW
      │   └── Mantidactylus_femoralis_CRH2340.fa *** NEW
      ├── /raw_data
      ├── File_rename.csv
      ├── /FrogCap_Files
      └── /Processed_Samples

*/ denotes directory

```

Once all the data is together in the "Assembled_Contigs" (or other folder), everything is ready for target marker matching. This step uses BLAST and the target marker fasta file to match each samples' contigs to the target markers. This script will detect and save paralogs separately and create a single large fasta file of all the samples' matching contigs renamed to the target marker for alignment in the second step. 

First, locate the target marker fasta file (i.e. "Master_Hyloidea_All-Markers_Apr21-2019.fa"), which will very soon all be publicly available. Second, edit the "03_Probe_Matching.R" R script, which will have this section: 

```r
##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#This script does the following:
#1. Matches the loci to the contigs, saves them to a new file
#2. Also finds the potential paralogs, removes them, and saves them to a separate file

#Set up directories
threads<-"6" #threads, keep quotes
contig.save<-"Project_Name"  #This is your save name for the big contig match file

#CLUSTER directories
work.dir<-"/home/username/Main_Project_Directory" #Your main project directory 
proc.dir<-"/home/username/Main_Project_Directory/Processed_Samples"
contig.dir<-"/home/username/Main_Project_Directory/Assembled_Contigs"
loci.file<-"/home/username/Main_Project_Directory/SELECT_YOUR_MARKER_FILE.fa"

```

When finished editing, run the script using the "Rscript" terminal command which installs alongside R. This can be done manually in the terminal window or added as the command in a cluster file. 

```bash
> Rscript 03_Probe_Matching.R

```

The marker matching script should take no longer than 3 hours unless you have a large number of samples. When it has finished, it should place two new files into each sample directory, in the "assembled-contigs" folder. The first, the "_orthologs.fa" file is the file of single-copy markers that will be used for alignment. The second file "_paralogs.fa" are potential paralogs, where more than 1 contig matched to the target marker. This could reflect poor assembly of that marker or a true paralog, and currently FrogCap does not use this file. In the future, FrogCap will include a paralog discovery tool to make these data usable. 

```
     /Spinomantis_elegans_CRH111
      ├── /reports
      ├── /assembled-contigs
      │   ├── Spinomantis_elegans_CRH111_contigs.fa
      │   ├── Spinomantis_elegans_CRH111_dipcontigs.fa
      │   ├── Spinomantis_elegans_CRH111_orthologs.fa *** NEW
      │   └── Spinomantis_elegans_CRH111_paralogs.fa *** NEW
      ├── /assembly-reads
      ├── /cleaned-reads-snp
      └── /raw-reads-symlink

*/ denotes directory

```

Once you confirm that these new files are present, return to the "Project_Name" directory to check that two new files are present. The first, the "-sample-assessment.csv" file provides summary statistics for each sample. These statistics include the number of matching orthologous markers, the number of paralogs, and matching contig length stats. The second file, the "Project_Name_contigs.fa" file is all the matching marker data extracted from each samples contigs. This fasta file includes the sequence data and the sample name and marker name in the header. 

```
     /Project_Name
      ├── /Assembled_Contigs
      ├── /raw_data
      ├── Project_Name_contigs.fa *** NEW
      ├── Project_Name-sample-assessment.csv *** NEW
      ├── File_rename.csv
      ├── /FrogCap_Files
      └── /Processed_Samples

*/ denotes directory

```

After running "03_Marker_Matching.R" and obtaining the new contigs file for alignment, move onto the next step to align each marker. 


## Target Marker Alignment [04_Marker_Alignment.R]

Insert explanation here for the next section 


```r
##########################################################################################################
#Parameter setups. Only edit values here. 
##########################################################################################################

#Parameters
threads<-"8"
min.taxa<-3 #min number taxa to keep an alignment

#Cluster dirs
species.loci<-"CONTIG-FILE-FROM-PREVIOUS-STEP.fa"  #The contig file output by the previous script 
work.dir<-"/home/username/Main_Project_Directory"
loci.file<-"/home/username/Main_Project_Directory/SELECT_YOUR_MARKER_FILE.fa" #Target loci file
out.dir<-"/home/username/Main_Project_Directory/Alignments" #The name of the output directory
```


When finished editing, run the script using the "Rscript" terminal command which installs alongside R. This can be done manually in the terminal window or added as the command in a cluster file. 

```bash
> Rscript 01_Pre_Process_Reads.R
```

Show example of output files and explain them. 

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

Once you confirm that the assembly finished properly for each sample, check the "Project_Name" directory to check that the "Assembled_Contigs" folder is present and contains all of the samples that were assembled in the previous step. This is also the point where additional contigs from other assembly runs can be added in, as the sample contigs placed in the "Assembled_Contigs" folder will be the samples used for sequence alignment. 

```
     /Project_Name
      ├── /Assembled_Contigs *** NEW
      │   ├── Spinomantis_elegans_CRH111.fa *** NEW
      │   ├── Boophis_burgeri_CRH0481.fa *** NEW
      │   ├── Aglyptodactylus_securifer_CRH1644.fa *** NEW
      │   └── Mantidactylus_femoralis_CRH2340.fa *** NEW
      ├── /raw_data
      ├── File_rename.csv
      ├── /FrogCap_Files
      └── /Processed_Samples

*/ denotes directory

```


When finished editing, run the script using the "Rscript" terminal command which installs alongside R. This can be done manually in the terminal window or added as the command in a cluster file. 

The Assembly run generated a statistics file in the "reports" folder and a new folder "assembled-contigs". This folder has two new files: the "_contigs.fa" file are the Spades assembled contigs while the "_dipcontigs.fa" file are the regular contigs where contigs from the same haplotype are merged together to create a "haplocontig". 

These two scripts can be run on any number of samples in any combination. However, when done with Tutorial 1, you will need to select your definitive sample set for alignment for Tutorial 2. 

"""

nbSave