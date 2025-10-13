# build the mutational model to calculate expected

python3 build_model/compile_denovo.py
python3 build_model/compare_denovo.py
python3 build_model/filter_denovo.py
python3 build_model/composite_likelihood_mito.py
python3 build_model/annotate_mutations.py

