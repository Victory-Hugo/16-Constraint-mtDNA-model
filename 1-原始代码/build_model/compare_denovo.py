from composite_likelihood_mito import make_type_count_vector, probability_per_ref_nuc, lr_ref_nuc, \
    probability_per_class, lr_class
import csv
import datetime
from typing import Dict, List, TextIO, Union

# global variables
nucleotides = ["A", "C", "G", "T"]
class_I_mutations = ["C>A", "T>A", "G>T", "A>T"]
class_II_mutations = ["C>G", "T>G", "G>C", "A>C"]
class_III_mutations = ["C>T", "T>C", "G>A", "A>G"]
mut_types = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]

# definitions of regions used in mitochondrial mutational model building
# ori refers to OriB-OriH region with known difference in mutational signature, m.16197-191, across artificial break
start_ori = 16197
end_ori = 191 + 1  # end not included in range so + 1
ori_region = list(range(start_ori, 16570)) + list(range(1, end_ori))
reference_except_ori = list(range(end_ori, start_ori))


def build_denovo_dict(variant: str, dict: Dict[str, int]):
    """Generate dictionary with the de novo variants and their total count.

    :param variant: de novo variant, in RefPosAlt format
    :param dict: dictionary name
    :return: dictionary where the variant is the key, and the value is the total count of the variant
    """
    if variant not in dict:
        dict[variant] = 1
    else:
        dict[variant] += 1
    return dict


def filter_denovo(germline_max: Union[int, None], som_tissue_max: Union[int, None], som_cancer_max: Union[int, None],
                  denovo_counts: Dict[str, int]):
    """Filters the all_denovo.txt file produced by compile_denovo.py.
    This is to remove de novo variants from samples that are above the maximum sample de novo count.

    :param germline_max: maximum sample de novo count for inclusion, for germline de novo
    :param som_tissue_max: maximum sample de novo count for inclusion, for somatic tissue de novo
    :param som_cancer_max: maximum sample de novo count for inclusion, for somatic cancer de novo
    :param denovo_counts: a dictionary where the variant is the key, and the value is the total count of the variant
    :return: denovo_counts dictionary
    """
    for row in csv.DictReader(open('output_files/denovo/all_denovo.txt'), delimiter='\t'):
        if "germline" in row["sample"]:
            if germline_max is not None:
                if int(row["sample_denovo_count"]) <= germline_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
        elif "somatic_tissue" in row["sample"]:
            if som_tissue_max is not None:
                if int(row["sample_denovo_count"]) <= som_tissue_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
        elif "somatic_cancer" in row["sample"]:
            if som_cancer_max is not None:
                if int(row["sample_denovo_count"]) <= som_cancer_max:
                    denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
            else:
                denovo_counts = build_denovo_dict(variant=row["denovo"], dict=denovo_counts)
    return denovo_counts


def compute_lr_classIII(denovo_counts: Dict[str, int], reference_region: List[int]):
    """Calculate the likelihood of mutation for the class III mutations.

    :param denovo_counts: a dictionary where the variant is the key, and the value is the total count of the variant
    :param reference_region: a list of the coordinates for the mtDNA region to analyze
    :return: lambda_classIII, a dictionary where the class III mutations are the keys and the value is their likelihood
    """
    mut_type_counts = make_type_count_vector(denovo_counts=denovo_counts, reference_region=reference_region)
    base_proportions = probability_per_ref_nuc(reference_region=reference_region)
    lambda_ref_nuc = lr_ref_nuc(type_counts=mut_type_counts, base_proportions=base_proportions)
    class_probabilities = probability_per_class(base_proportions=base_proportions)
    lambda_mut_class = lr_class(type_counts=mut_type_counts, class_probabilities=class_probabilities)

    lambda_classIII = {}
    for mut in class_III_mutations:
        lambda_classIII[mut] = lambda_ref_nuc[mut[0]] * lambda_mut_class[mut]
    return lambda_classIII


def apply_threshold(germline_max: int, som_tissue_max: int, som_cancer_max: int, category: str, threshold: int,
                    reference_region: List[int], file: TextIO):
    """For a category of de novo, iterate through a range of thresholds for the maximum sample de novo count.
    Any de novo variants from samples with a de novo count greater than threshold will not be included in calculations
    for the class III likelihoods.

    :param germline_max: maximum sample de novo count for inclusion, for germline de novo
    :param som_tissue_max: maximum sample de novo count for inclusion, for somatic tissue de novo
    :param som_cancer_max: maximum sample de novo count for inclusion, for somatic cancer de novo
    :param category: the category of de novo where the maximum sample de novo count is being iteratively increased
    :param threshold: the maximum sample de novo mutations count in the loop, for print to file
    :param reference_region: a list of the coordinates for the mtDNA region to analyze
    :param file: the txt file to write the class III likelihoods to
    """
    denovo_counts = {}
    denovo_counts = filter_denovo(germline_max=germline_max, som_tissue_max=som_tissue_max,
                                  som_cancer_max=som_cancer_max, denovo_counts=denovo_counts)

    print('\n' + "For ", category, " testing maximum sample de novo count of ", threshold, '\n')
    lambda_classIII = compute_lr_classIII(denovo_counts=denovo_counts, reference_region=reference_region)

    file.write(str(category) + '\t' + str(threshold) + '\t' +
               str(lambda_classIII["C>T"]) + '\t' + str(lambda_classIII["G>A"]) + '\t' +
               str(lambda_classIII["T>C"]) + '\t' + str(lambda_classIII["A>G"]) + '\n')


if __name__ == "__main__":
    print(datetime.datetime.now(), "Starting to compare de novo categories!")
    print('\n' + "This will produce likelihood scores for transitions across sample categories")

    f = open('output_files/denovo/Ts_likelihood_by_category.txt', "w")
    f.write("category	threshold	likelihood_C>T	likelihood_G>A	likelihood_T>C	likelihood_A>G" + '\n')

    # range of values to use as the maximum sample de novo mutations count, for inclusion
    threshold_range = [1, 2, 3, 4, 5, 40]  # I know maximum is 38, so have manually determined the thresholds to test

    for threshold in threshold_range:
        # first, iterate through germline only
        apply_threshold(germline_max=threshold, som_tissue_max=0, som_cancer_max=0,
                        category="germline", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # then, include all germline plus iterate through somatic tissue
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=threshold, som_cancer_max=0,
                        category="germline + somatic tissue", threshold=threshold,
                        reference_region=reference_except_ori, file=f)

        # then, include all germline and somatic tissue, plus iterate through somatic cancer
        apply_threshold(germline_max=max(threshold_range), som_tissue_max=max(threshold_range),
                        som_cancer_max=threshold, category="germline + somatic tissue + somatic cancer",
                        threshold=threshold, reference_region=reference_except_ori, file=f)

    print(datetime.datetime.now(), "Script finished!")
