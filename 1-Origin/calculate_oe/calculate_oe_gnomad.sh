# calculate the oe ratio for functional classes of variation across the mtDNA and within genes in gnomAD
# these scripts have been set up such that the parameters for gnomAD analyses are defaults and don't need to be provided

python3 calculate_oe/calibrate_model_step1.py
Rscript calculate_oe/calibrate_model_step2.R
python3 calculate_oe/oe_consequences.py
python3 calculate_oe/oe_loci.py
python3 calculate_oe/oe_other.py
