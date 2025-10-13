## Overview:

These scripts were used to build the novel mitochondrial mutational model.

`compile_denovo.py`: Parse lists of de novo variants obtained from the literature and an in-house dataset. Outputs list of all de novo calls and their source. The datasets ascertained from publications are provided in `required_files`; **note that the one unpublished dataset used is not provided in this repo due to data restrictions, thus results will differ slightly when only using the available datasets in this repo**. Details on the de novo datasets used is provided in the manuscript Supplementary Information. 

`compare_denovo.py`: Compare the mutational likelihoods of transitions, across different sources (germline, somatic tissue, somatic cancer) and sample de novo counts.

`filter_denovo.py`: Remove any de novo from outlier samples. Outputs list of de novo mutations used to calculate mutability. User specified arguments *-germline_max*, *-som_tissue_max*, *-som_cancer_max* to indicate the maximum sample de novo count used for filtering for each source. Defaults provided based on analyses.

`composite_likelihood_mito.py`: Code to apply the mitochondrial composite likelihood model. Outputs file with mutational likelihood scores for every single nucleotide variant in the mtDNA. User specified arguments *-context_size* to indicate how many nucleotides to include for sequence context, 3 for trinucleotide set as default, and *-denovo_list* the path to the list of de novo to use, the output of `filter_denovo.py` set as default.

`annotate_mutations.py`: Annotate the output of `composite_likelihood_mito.py` with annotations needed for downstream analyses, such as variant consequence, gene, and in silico predictions. User specified argument *-input*; the output of `composite_likelihood_mito.py` set as default.

`build_model.sh`: Apply the above analyses to produce a list of annotated mtDNA variants and their likelihood scores. Run using````bash -e build_model/build_model.sh````.

`simulate_heteroplasmy.R`: Apply a computational model of germline mtDNA mutation and heteroplasmy drift to support a correlation between mutation rates and maximum heteroplasmy.

#### Example running simulation:

Apply a computational model of germline mtDNA mutation and heteroplasmy drift. Run time <5 minutes using 12 cores.

````
Rscript simulate_heteroplasmy.R
````
Produces text files outputs that are used for `figure_scripts/FigureS4.R`.