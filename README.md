# FrogCap Sequence Capture probeset and pipeline

Pipeline for the analysis of data acquired from probe set designed by Hutter et al. (2019) in review and available as a preprint:

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and Brown R.M. (2019). FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. bioRxiv 825307. doi: https://doi.org/10.1101/825307


The pipeline here is in the process of being upgraded to version 2, and utilizes many open source programs to analyze genomic sequence capture data. 

The probe set files themselves will be publicly available once the manuscript for the probe set is published. However, if you would like to collaborate, the following probe sets are available: 

1) Amphibian-level probe set: this probe set will function across any group of amphibians and maintains compatibility across Amphibia. The disadvantage is that variability is lost due to the conserved nature of the markers, so a more narrow probe set might be best for smaller phylogenetic questions. This probe set is also best for Archaeobatrachian (i.e. basal) frog lineages (Not yet published, available upon request). 


2) Ranoidea superfamily probe-set: This probe set will function across the Ranoidea superfamily of frogs. While some ultra-conserved exons and elements are included that are compatible across frogs, the main priority of this probe set is to acquire exons that are variable within Ranoidea so that more fine-scale phylogenetic questions can be addressed. 

		Ranoidea V1: the probe set used in Hutter et al. 2019

 		Ranoidea V2: removes markers not represented well throughout the superfamily. 
		Adds in more UCEs, Legacy markers. Recommended. 


3) Hyloidea superfamily probe-set: This probe set will function across the Hyloidea superfamily of frogs. While some ultra-conserved exons and elements are included that are compatible across frogs, the main priority of this probe set is to acquire exons that are variable within Hyloidea so that more fine-scale phylogenetic questions can be addressed. 

		Hyloidea V1: the probe set used in Hutter et al. 2019

		Hyloidea V2: removes markers not represented well throughout the superfamily. 
		Adds in additional shared markers with Ranoidea. Recommended. 


All three probe sets contain 2200 UCEs from Streicher et al. 2018. In addition, about ~6000 exons are shared across the probe sets, with different probes that will work more effectively depending on the group chosen above. Additionally, custom probe sets for your target group of interest can be created upon request, and a web form will soon be available to run the customization yourself.

Please contact me if you have any questions! 

## Tutorials 

The Wiki above has instructions and tutorials for using all the scripts here! Quick links:

[Tutorial 0: Customized FrogCap probeset](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-0:-Customized-FrogCap-probe-set)

[Tutorial 1: Cluster and program setup](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-1---Setting-up-environment)

[Tutorial 2: Data preprocessing (cleaning raw reads, decontamination, assembly)](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-2---Data-Processing)

[Tutorial 3: Data alignment (target marker matching, alignment, alignment filtering)](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Tutorial-3---Marker-Matching-and-Alignment)

Tutorial 4: Tree estimation (gene trees, concatenation, species trees) COMING SOON

[Troubleshooting: random issues that arise](https://github.com/chutter/FrogCap-Sequence-Capture/wiki/Troubleshooting)


## Papers that have used FrogCap:

Rasolonjatovo S.M., Scherz M.D., Hutter C.R., Glaw F., Rakotoarison A., Razafindraibe J.H., Goodman S.M., Raselimanana A.P. and Vences M. (2020). Sympatric lineages in the Mantidactylus ambreensis complex of Malagasy frogs originated allopatrically rather than by in-situ speciation. Molecular Phylogenetics and Evolution, 144, 106700.

Hutter C.R., Cobb K.A., Portik D., Travers S., Wood Jr. P.L., and Brown R.M. (2019). FrogCap: A modular sequence capture probe set for phylogenomics and population genetics for Anurans, assessed across multiple phylogenetic scales. bioRxiv 825307. doi: https://doi.org/10.1101/825307

Chan K.O., Hutter C.R., Wood P.L., Grismer L.L., and Brown R.M. (2019). Exons, Introns, and UCEs Reveal Conflicting Phylogenomic Signals in a Rapid Radiation of Frogs (Ranidae: Hylarana). bioRxiv, p.765610. doi: https://doi.org/10.1101/765610

Chan K.O., Hutter C.R., Wood P.L., Grismer L.L., and Brown R.M. (2019). Species delimitation in the grey zone: introgression obfuscates phylogenetic inference and species boundaries in a cryptic frog complex (Ranidae: Pulchrana picturata). bioRxiv, 832683. doi: https://doi.org/10.1101/832683



