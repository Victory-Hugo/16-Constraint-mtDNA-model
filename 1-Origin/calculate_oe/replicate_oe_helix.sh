# calculate the oe ratio for functional classes of variation in a replication dataset with heteroplasmy, HelixMTdb
# HelixMTdb excludes variants at these bases: 300-316, 513-525, and 16182-16194
# helix_excluded = list(range(300, 317)) + list(range(513, 526)) + list(range(16182, 16195))

mkdir output_files/oe/replication_dataset/
mkdir output_files/calibration/replication_dataset/
python3 calculate_oe/calibrate_model_step1.py -obs 'helix_max_hl' -prefix 'replication_dataset/helix_' -exc_sites 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 513 514 515 516 517 518 519 520 521 522 523 524 525 16182 16183 16184 16185 16186 16187 16188 16189 16190 16191 16192 16193 16194
Rscript calculate_oe/calibrate_model_step2.R 'replication_dataset/helix_' 'no'
python3 calculate_oe/oe_consequences.py -obs 'helix_max_hl' -parameters 'output_files/calibration/replication_dataset/helix_linear_model_fits.txt' -prefix 'replication_dataset/helix_' -exc_sites 300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 513 514 515 516 517 518 519 520 521 522 523 524 525 16182 16183 16184 16185 16186 16187 16188 16189 16190 16191 16192 16193 16194


