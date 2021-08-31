# FrogCap Sequence Capture probe set and pipeline

Pipeline for the analysis of data acquired from probe set designed by Hutter et al. (2019; in review) and available as a preprint:

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and Brown R.M. (2019). FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. bioRxiv 825307. doi: https://doi.org/10.1101/825307

The data analysis pipeline above (scripts 1-4) utilizes many open source programs to analyze genomic sequence capture data and is described in Hutter et al. (2019). The scripts need very little technical knowledge, and require some parameters to be set. Below there are tutorials for setting up and executing each step of the pipeline process.  

The probe set files are publicly available above in the Probe Sets folder. The probe sets are:

1) Ranoidea superfamily probe set: This probe set will function across the Ranoidea superfamily of frogs. While some ultra-conserved exons and elements are included that are compatible across frogs, the main priority of this probe set is to acquire exons that are variable within Ranoidea so that more fine-scale phylogenetic questions can be addressed. 

		Ranoidea V1: the probe set used in Hutter et al. (2019)

 		Ranoidea V2: removes markers not represented well throughout the superfamily. 
		Adds in more UCEs, Legacy markers. Recommended. 


2) Hyloidea superfamily probe set: This probe set will function across the Hyloidea superfamily of frogs. While some ultra-conserved exons and elements are included that are compatible across frogs, the main priority of this probe set is to acquire exons that are variable within Hyloidea so that more fine-scale phylogenetic questions can be addressed. 

		Hyloidea V1: the probe set used in Hutter et al. (2019)

		Hyloidea V2: removes markers not represented well throughout the superfamily. 
		Adds in additional shared markers with Ranoidea. Recommended. 

3) Amphibian-level probe set: this probe set will function across any group of amphibians and maintains compatibility across Amphibia. The disadvantage is that variability is lost due to the conserved nature of the markers, so a more narrow probe set might be best for smaller phylogenetic questions. This probe set is also best for Archaeobatrachian (i.e. basal) frog lineages. This probe set also contains some BUSCO markers to allow compatibility with other vertebrate probe sets. This probe set is the only one not yet published, collaboration is available upon request. 

All three probe sets contain 2200 UCEs from Streicher et al. 2018. In addition, about ~6000 exons are shared across the probe sets, with different probes that will work more effectively depending on the group chosen above. Additionally, custom probe sets for your target group of interest can be created upon request, or you may use the script in Tutorial 0 to design a custom probe set. 

Please contact me if you have any questions! 

## Tutorials 

The Wiki above has instructions and tutorials for using all the scripts here! Quick links:

[Tutorial 0: Customized FrogCap probeset](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-0:-Customized-FrogCap-probe-set)

[Tutorial 1: Cluster and program setup](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-1---Setting-up-environment)

[Tutorial 2: Data preprocessing (cleaning raw reads, decontamination, assembly)](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-2---Data-Processing)

[Tutorial 3: Data alignment (target marker matching, alignment, alignment filtering)](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-3---Marker-Matching-and-Alignment)

[Troubleshooting: random issues that arise](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Troubleshooting)


## Papers that have used FrogCap:

Rasolonjatovo S.M., Scherz M.D., Hutter C.R., Glaw F., Rakotoarison A., Razafindraibe J.H., Goodman S.M., Raselimanana A.P. and Vences M. (2020). Sympatric lineages in the Mantidactylus ambreensis complex of Malagasy frogs originated allopatrically rather than by in-situ speciation. Molecular Phylogenetics and Evolution, 144, 106700.

Chan K.O., Hutter C.R., Wood P.L., Grismer L.L., and Brown R.M. (2019). Exons, Introns, and UCEs Reveal Conflicting Phylogenomic Signals in a Rapid Radiation of Frogs (Ranidae: Hylarana). bioRxiv, p.765610. doi: https://doi.org/10.1101/765610

Chan K.O., Hutter C.R., Wood P.L., Grismer L.L., and Brown R.M. (2019). Species delimitation in the grey zone: introgression obfuscates phylogenetic inference and species boundaries in a cryptic frog complex (Ranidae: Pulchrana picturata). bioRxiv, 832683. doi: https://doi.org/10.1101/832683

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and Brown R.M. (2019). FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. bioRxiv 825307. doi: https://doi.org/10.1101/825307
