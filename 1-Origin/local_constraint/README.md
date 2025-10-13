## Overview:

`local_constraint.py`: Calculate Mitochondrial Local Constraint (MLC) scores. User specified argument *-input*;  the output of `annotate_mutations.py` set as default, *-kmer_length*; 30 base pairs set as default. All other arguments set as defaults for gnomAD (*-obs*, *-parameters*, *-prefix*, *-exc_sites*). Outputs *kmers_local_constraint.txt* with annotations for every kmer, and *per_base_local_constraint.txt* with annotations for every base.

Uses multiprocessing module to maximise use of available CPU.