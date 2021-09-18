import nimib, nimibook

nbInit

nbUseNimibook 

nbText: md"""
# Create your own customized FrogCap probe set

This page will guide you on creating your own customized FrogCap probe set. You can draw markers from the growing FrogCap database that have been previously tested, and thus the output baits are ready for synthesis with no further processing necessary. Functionality includes: (1) The ability to select markers from certain families, subfamilies, and genera; (2) Different types of markers (UCEs, exons); (3) and marker features such as proportion of informative sites and marker length. 

To begin, download the "Custom-Probe-Set" folder in the main GitHub directory to your computer. The folder should contain the following files:

```
     /Custom-Probe-Set
      ├── FrogCap_Probe-Select.R
      ├── Taxonomy_dataset.csv
      └── Marker_Bait_dataset.csv

 *** newly generated files
```

The folder contains the three files:

File | Explanation
-----------|:--------:|
FrogCap_Probe-Select.R | The R script for creating the custom probe set
Taxonomy_dataset.csv | A taxonomy database for each marker and taxonomic group sequenced
Marker_Bait_dataset.csv | Markers and their associated baits and other metadata


Next, open the folder and edit the R script titled: "FrogCap_Probe-Select.R", which will appear like so:


```r
###############################################################################
########   Step 1: Create locus summary dataset               #################
###############################################################################
###############################################################################
###############################################################################

###### Directory setup ########

work.dir = "/Full/Directory/Path/to/Custom-Probe-Set"
output.file ="Custom"

#The two metadata files should be placed in the working directory 
tax.data.file = "Taxonomy_dataset.csv"
marker.data.file = "Marker_Bait_dataset.csv"

#################################################
#Features you would like the probe set to have
#################################################

probe.size = 20000 #Number of probes, maximum is currently 20k for custom kits. use full kit
uce.include = TRUE #Other option FALSE, no UCEs
most.informative = TRUE #FALSE = randomly selects markers rather than by information content
min.exon.length = 200 #The minimum exon length
max.exon.length = 100000 #The maximum exon length

#Select a family or subfamily from the two lists below
clade.select = "Microhylidae"

#You can also select two or more taxonomic groups [overrides previous clade select]
#The program will select markers present in BOTH. 
#Must both be the same phylogenetic scale (Family / Subfamily)
#Cannot be from different Superfamilies (Ranoidea / Hyloidea)
clade.select = c("Mantellidae", "Rhacophoridae")

```

This section is the only section that needs editing. Namely, the directories should be the full path to make things go smoother. You can set various parameters here for the probe set. When finished editing, run the R script in the command line:

```bash

> Rscript FrogCap_Probe-Select.R

```

Which should have generated 3 files: 


```
     /Custom-Probe-Set
      ├── FrogCap_Probe-Select.R
      ├── Taxonomy_dataset.csv
      ├── Marker_Bait_dataset.csv
      ├── Custom_Baits.fa ***
      ├── Custom_Markers.fa ***
      └── Custom_metadata.txt ***

 *** newly generated files
```

And the three output files are: (1) Custom_baits.fa, which are the bait sequences as a separate Fasta entry per bait. These are ready for synthesis and require no further processing since this has already been done; (2) Custom_Markers.fa, which are the assembled markers the baits are derived from. This is used in downstream pipeline analyses; and (3) Custom_metadata.txt, which are the metadata that defined each marker for your information. Not needed for any analysis. 

"""

nbSave