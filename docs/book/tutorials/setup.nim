import nimib, nimibook

nbInit

nbUseNimibook 

nbText: md"""
# Setting up the environment

## Required software 

The FrogCap pipeline uses R as a wrapper to run multiple programs across all of the samples in the dataset, rather than one at a time. It also organizes and moves files around, deletes temporary files, and keeps your data and workflow organized. Later tutorials in the pipeline will use specific packages from R and BIOCONDUCTOR. 

***

1. Essentially a program and package manager. Especially useful for switching between different environments where programs can have conflicting dependencies. Also makes using R on a cluster easier.

* **anaconda**: https://www.anaconda.com/distribution/

2. A program that matches sequence data to other sequence data.

* **BLAST**: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Everything below can be installed via Anaconda. Links are provided in case the Anaconda install does not work or Anaconda is not being used. Currently (May 10 2019) AMAS and BLAST should be installed manually as they do not function properly through Anaconda.


3. Removal of adapter contamination and low quality bases:

* **Fastp**: https://github.com/OpenGene/fastp

4. Set of tools including short read matching and other processing tools: 

* **bbtools**: https://jgi.doe.gov/data-and-tools/bbtools/

5. Program that assembles short-read Illumina read data into contigs: 

* **Spades**: http://cab.spbu.ru/software/spades/

6. Alignment of sequence data: 

* **MAFFT V7**: https://mafft.cbrc.jp/alignment/software/

7. Concatenation and data manipulation: 

* **AMAS**: https://github.com/marekborowiec/AMAS

8. Trimming of sequence alignments: 

* **TRIMAL**: http://trimal.cgenomics.org

9. The main scripting language used for this pipeline: 

* **R statistical software**: https://www.r-project.org

Required core R packages
* ape
* seqinr
* stringr
* data.table

Required BIOCONDUCTOR R packages
* GenomicRanges
* Biostrings
* Rsamtools


## Download and organize required files

DOWNLOAD LINK COMING SOON (available upon request also)

Note the probe set and version you are using in the "Probe_Sets" folder, and use those locus files only. The "Contamination_Genomes" folder contains genomes of potential contaminants from outside organisms. If you do not wish to use these, just delete or not include the folder and the scripts will ignore that step (not recommended).  


## Create new Anaconda environment

On clusters, an efficient method of managing software is through Anaconda. Anaconda allows separate environments to be created for specific purposes, as dependencies can conflict and cause everything to stop working. 

After installing Anaconda, create a new environment as shown below. "--prefix" is the path to your environment, where here I chose $WORK because Anaconda generates a large number of files and $HOME has a file number limit. If your cluster has the same restriction, you should also install everything in $WORK (or elsewhere). "conda" is next the folder that will contain several different environments, and "FrogCap" will be the name of this current environment. Finally, "source activate" will activate the new environment, so that new programs can be installed. 

>
> conda create --prefix $WORK/conda/FrogCap
>
> source activate $WORK/conda/FrogCap
>

Next, set up the channels where Anaconda will look for programs to be downloaded. Input the "conda config" commands in that order, as multiple programs are available in different repositories, but can sometimes be outdated. This order ensures the most recent versions are installed. 

>
> conda config --add channels bioconda
>
> conda config --add channels conda-forge
>
> conda config --add channels defaults
>

Install all the programs available on Anaconda as shown here: 

>
> conda install -c r r-base
>
> conda install -c bioconda amas fastp mafft trimal bbmap
>

Finally, install all the R packages required by FrogCap: 

>
> conda install -c r r-stringr
>
> conda install -c conda-forge r-ape r-seqinr r-data.table
>
> conda install -c bioconda bioconductor-genomicranges bioconductor-biostrings bioconductor-rsamtools
>

Note that programs should be grouped together as shown to ensure that the correct dependencies are installed and not later overwritten. 


## Download and install programs to your $PATH 

Note: Needed for Blast and AMAS since their Anaconda installations do not function properly. 

Some programs may need to be downloaded and manually compiled and placed into your $PATH variable if Anaconda fails to install it or a dependency is changed and breaks the program. Alternatively, if Anaconda is not being used, you will have to do this manually for each program. On OSX and Linux, there should be a file called ".bashrc" or ".bash_profile" that has a list of your programs and their absolute paths (i.e. location on the computer), which the terminal uses to locate the program when its name is entered. 

For example, in my $HOME directory on our cluster, I have a folder called "programs" to store programs I have downloaded not available on anaconda. 

In particular, you should manually instal Blast, which for some reason breaks R in anaconda:

After downloading and unzipping the Blast download folder, open the ".bashrc" or ".bash_profile" file and enter:

>
> PATH="/home/c111h652/programs/blast-path:$PATH"
>

Where the "blast-path" folder contains all the different Blast executables. After restarting your terminal session (the bash file loads on terminal startup), you should now be able to run Blast from any directory. 

"""

nbSave