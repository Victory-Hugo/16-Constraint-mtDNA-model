## Overview:

`calibrate_model_step1.py`: For each gene/locus, calculate the observed sum maximum heteroplasmy of neutral variants and their sum mutation likelihood scores to fit linear models. User specified argument *-input*, the file with likelihood scores and observed values; the output of `annotate_mutations.py` set as default.

`calibrate_model_step2.R`: Use the output of `calibrate_model_step1.py` to fit linear models for netural variation. 

`oe_functions.py`: Functions used for observed:expected value calculations.

`oe_consequences.py`: Calculate the observed:expected ratios and confidence intervals for different functional classes of variation. User specified argument *-input*;  the output of `annotate_mutations.py` set as default. All other arguments set as defaults for gnomAD (*-obs*, *-parameters*, *-prefix*, *-exc_sites*). 

`oe_loci.py`: Calculate the observed:expected ratios and confidence intervals for different functional classes of variation in different genes/loci. User specified argument *-input*;  the output of `annotate_mutations.py` set as default. All other arguments set as defaults for gnomAD (*-obs*, *-parameters*, *-prefix*, *-exc_sites*). 

`oe_other.py`: Calculate the observed:expected ratios and confidence intervals for other categories of variation. User specified argument *-input*;  the output of `annotate_mutations.py` set as default. All other arguments set as defaults for gnomAD (*-obs*, *-parameters*, *-prefix*, *-exc_sites*). 

`calculate_oe_gnomad.sh`: Apply the above analyses to assess constraint in gnomAD.

`replicate_oe_helix.sh`: Apply the above analyses to assess constraint in a replication dataset, HelixMTdb. User specified arguments provided for HelixMTdb to override default of gnomAD.

`oe_lookup.py`: Calculate the observed:expected ratio and confidence interval between a pair of coordinates, using gnomAD.

