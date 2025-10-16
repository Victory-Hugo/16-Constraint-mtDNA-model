#!/bin/bash
#!不建议一次性运行，因为代码中有绝对路径
#!千万不要修改python代码文件名
#!如下路径千万不要修改，也不要重命名文件名、文件路径等:



cd /mnt/f/OneDrive/文档（科研）/脚本/Download/16-Constraint-mtDNA-model/2-群体版本

python3 2-计算观测值和预期值/calibrate_model_step1.py
Rscript 2-计算观测值和预期值/calibrate_model_step2.R
python3 2-计算观测值和预期值/oe_consequences.py
python3 2-计算观测值和预期值/oe_loci.py
python3 2-计算观测值和预期值/oe_other.py

#* 如果需要计算,请给出起始和终止位置：
python3 \
    2-计算观测值和预期值/oe_lookup.py \
    -start 1000 \
    -end 2000




mkdir -p output/oe/replication_dataset/
mkdir -p output/calibration/replication_dataset/


python3 2-计算观测值和预期值/calibrate_model_step1.py \
    -obs 'helix_max_hl' \
    -prefix 'replication_dataset/helix_' \
    -exc_sites 300 301 302 303 304 305 306 307 308 309 310 \
    311 312 313 314 315 316 513 514 515 516 517 518 519 520 \
    521 522 523 524 525 16182 16183 16184 16185 16186 16187 \
    16188 16189 16190 16191 16192 16193 16194

Rscript 2-计算观测值和预期值/calibrate_model_step2.R \
    'replication_dataset/helix_' \
    'no'

python3 2-计算观测值和预期值/oe_consequences.py \
    -obs 'helix_max_hl' \
    -parameters 'output/calibration/replication_dataset/helix_linear_model_fits.txt' \
    -prefix 'replication_dataset/helix_' \
    -exc_sites 300 301 302 303 304 305 306 307 308 309 310 \
    311 312 313 314 315 316 513 514 515 516 517 518 519 520 \
    521 522 523 524 525 16182 16183 16184 16185 16186 16187 \
    16188 16189 16190 16191 16192 16193 16194
