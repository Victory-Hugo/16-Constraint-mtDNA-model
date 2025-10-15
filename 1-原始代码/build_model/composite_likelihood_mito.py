import argparse
from compile_denovo import rcrs_pos_to_ref
import csv
import datetime
from typing import Dict, List, Tuple, Union
import os

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


# helper functions
def build_dictionary(key: Union[str, Tuple[str, str, int]], dictionary: Dict[str, int]):
    """Generate a dictionary that provides a count for the key.

    :param key: the string or tuple to count occurrences of
    :param dictionary: the name of the dictionary
    :return: a dictionary, where the key is the string or tuple to count occurrences of, and the value is the count
    """
    if key not in dictionary:
        dictionary[key] = 1
    else:
        dictionary[key] += 1
    return dictionary


def build_additive_dictionary(key: Union[str, Tuple[str, str, int]], value: int, dictionary: Dict[str, int]):
    """Generate a dictionary that provides a cumulative count for the key, where each key has a non-zero value.

    :param key: the string or tuple to sum the value across
    :param value: the count or value assigned to the key
    :param dictionary: the name of the dictionary
    :return: a dictionary, where the key is the string or tuple to sum across, and the value is the cumulative count
    """
    if key not in dictionary:
        dictionary[key] = value
    else:
        dictionary[key] += value
    return dictionary


def flip_base(base: str):
    """Return the pairing nucleotide that is located in the reverse complement strand at the same position.

    :param base: the reference nucleotide
    :return: flipped_base, the pairing nucleotide
    """
    flipped_base = ''
    if base == "A":
        flipped_base = "T"
    elif base == "C":
        flipped_base = "G"
    elif base == "G":
        flipped_base = "C"
    elif base == "T":
        flipped_base = "A"
    return flipped_base


def flip_mut(mut: str):
    """Return the equivalent mutation for the nucleotide at the same position on the reverse complement strand.

    :param mut: the mutation type in Ref>Alt format
    :return: flipped_mut, the equivalent mutation for the nucleotide on the reverse complement strand
    """
    flipped_mut = ''
    if mut == "A>C":
        flipped_mut = "T>G"
    elif mut == "A>G":
        flipped_mut = "T>C"
    elif mut == "A>T":
        flipped_mut = "T>A"
    elif mut == "G>A":
        flipped_mut = "C>T"
    elif mut == "G>C":
        flipped_mut = "C>G"
    elif mut == "G>T":
        flipped_mut = "C>A"
    elif mut == "C>A":
        flipped_mut = "G>T"
    elif mut == "C>G":
        flipped_mut = "G>C"
    elif mut == "C>T":
        flipped_mut = "G>A"
    elif mut == "T>A":
        flipped_mut = "A>T"
    elif mut == "T>C":
        flipped_mut = "A>G"
    elif mut == "T>G":
        flipped_mut = "A>C"
    return flipped_mut


def handle_mito_genome(flanking_coord: int, coord: int):
    """For flanking nucleotides around the artificial break in the circular mtDNA, or around m.3107N spacer, return
    the correct flanking nucleotide coordinate to use. Can handle correction of flanking nucleotide range [-10:10]/{0}.

    :param flanking_coord: the coordinate of the flanking nucleotide, before correction
    :param coord: the coordinate of the reference nucleotide
    :return: the correct coordinate of the flanking nucleotide as string, note this will be unchanged unless corrected
    """
    if flanking_coord < 1:  # to handle circular genome
        flanking_coord = 16569 + flanking_coord  # since flanking_coord will be 0 or negative number
    elif flanking_coord > 16569:  # to handle circular genome
        flanking_coord = flanking_coord - 16569
    elif flanking_coord < 3107 and 3118 > coord > 3107:  # to handle m.3107N, this skips this position
        flanking_coord = flanking_coord - 1
    elif flanking_coord > 3107 and 3096 < coord < 3107:  # to handle m.3107N, this skips this position
        flanking_coord = flanking_coord + 1
    elif flanking_coord == 3107 and flanking_coord < coord:  # to handle m.3107N
        flanking_coord = 3106
    elif flanking_coord == 3107 and flanking_coord > coord:  # to handle m.3107N
        flanking_coord = 3108
    flanking_coord = str(flanking_coord)
    return flanking_coord


def write_file_header(file: str, flanking_range: List[int], variable: str):
    """Creates the file and header to write sequence context related data to file, for plotting.

    :param file: the path of the file to write to
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param variable: the data point to compare the flanking nucleotide to
    :return: f, the file to which the data is being written
    """
    f = open(file, "w")
    f.write("flanking_nucleotide" + '\t')
    for flanking_pos in flanking_range:
        if flanking_pos == flanking_range[-1]:
            f.write(str(flanking_pos) + '\t' + variable + '\n')
        else:
            f.write(str(flanking_pos) + '\t')
    return f


# functions for the composite likelihood model
def make_denovo_counts(denovo_list: str):
    """Read in the final de novo variants to use, and convert to dictionary.

    :param denovo_list: path to the final list of de novo mutations to use, in RefPosAlt format with their counts
    :return: denovo_counts, a dictionary where the variant is the key, and the value is the total count of the variant
    """
    denovo_counts = {}
    for row in csv.DictReader(open(denovo_list), delimiter='\t'):
        denovo_counts[row["denovo"]] = int(row["count"])
    return denovo_counts


def make_type_count_vector(denovo_counts: Dict[str, int], reference_region: List[int]):
    """Count the number of each mutation type within the region specified.

    :param denovo_counts: a dictionary where the variant is the key, and the value is the total count of the variant
    :param reference_region: a list of coordinates for the mtDNA region to count the occurrence of each mutation within
    :return: mut_type_counts, dictionary where key is the mutation type and the value is its total count in the region
    """
    mut_type_counts = {}
    # use a list of all possible SNVs to parse quickly
    for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        variant = row["REF"] + row["POS"] + row["ALT"]
        mut_type = row["REF"] + ">" + row["ALT"]
        coord = row["POS"]
        # count the number of each 12 types of mutation
        if variant in denovo_counts:
            if int(coord) in reference_region:
                allele_count = denovo_counts[variant]
                mut_type_counts = build_additive_dictionary(key=mut_type, value=allele_count,
                                                            dictionary=mut_type_counts)
    print('\n' + "Number of observed de novo mutations: ", sum(mut_type_counts.values()))
    print('\n' + "Mutation type counts for region: ", mut_type_counts)
    return mut_type_counts


def probability_per_ref_nuc(reference_region: List[int]):
    """Count the number of times each base (A,C,G,T) occurs within a mtDNA region, and express as a proportion.
    This will be used to represent the probability of mutation at each reference nucleotide.

    :param reference_region: a list of the coordinates for the mtDNA region to count the occurrence of each base within
    :return: base_proportions, a dictionary where key is the base and the value is the proportion of base in the region
    """
    # use fasta file of the reference mitochondrial genome sequence, use replace to enable parsing per position
    fasta = open('required_files/synthetic_vcf/NC_012920.1_noheader.fasta').read().replace('\n', '')

    dictionary = {}
    # build a dictionary to count the number of A, C, G, T in the region
    for coord in reference_region:
        # require the -1 to convert coordinate to position in fasta since the first position in the fasta is 0, not 1
        if fasta[coord - 1] != "N":  # if base is not N, to skip position m.3107
            dictionary = build_dictionary(key=fasta[coord - 1], dictionary=dictionary)
    total_number_bases = sum(dictionary.values())

    # calculate as a proportion, keys = A, C, G, T
    base_proportions = {}
    for key in dictionary:
        base_proportions[key] = dictionary[key] / total_number_bases
    print('\n' + "The proportion of each base is: ", base_proportions)
    return base_proportions


def lr_ref_nuc(type_counts: Dict[str, int], base_proportions: Dict[str, float]):
    """Calculate the likelihood ratio of mutation at each reference nucleotide n(t) = A,C,G,T.

    :param type_counts: dictionary where key is the mutation type and value is the total count of type
    :param base_proportions: a dictionary where key is the base and the value is the proportion of the base
    :return: lambda_ref_nuc, where key is the reference nucleotide and value is likelihood ratio of reference nucleotide
    """
    sum_vtype = sum(type_counts.values())
    lambda_ref_nuc = {}
    # for each base, tally the number of mutations where it is the reference nucleotide, then calculate likelihood ratio
    for base in nucleotides:
        total_mut_at_base = 0
        for key in type_counts:
            if key.startswith(base + ">"):
                total_mut_at_base += type_counts[key]
        lambda_ref_nuc[base] = total_mut_at_base / (sum_vtype * base_proportions[base])
    print('\n' + "The likelihood of mutation at each reference nucleotide is: ", lambda_ref_nuc)
    return lambda_ref_nuc


def probability_per_class(base_proportions: Dict[str, float]):
    """Calculate the probability of each mutation class, using information about the frequency of alternate nucleotide.

    :param base_proportions: a dictionary where key is the base and the value is the proportion of base
    :return: class_probabilities, dictionary where key is class and value is the probability of the mutation class
    """
    class_probabilities = {}
    probC_to_A = base_proportions["C"] * (base_proportions["A"] / (1 - base_proportions["C"]))  # type 1
    probG_to_T = base_proportions["G"] * (base_proportions["T"] / (1 - base_proportions["G"]))  # type 7
    probT_to_A = base_proportions["T"] * (base_proportions["A"] / (1 - base_proportions["T"]))  # type 4
    probA_to_T = base_proportions["A"] * (base_proportions["T"] / (1 - base_proportions["A"]))  # type 10
    class_probabilities["I"] = probC_to_A + probG_to_T + probT_to_A + probA_to_T

    probC_to_G = base_proportions["C"] * (base_proportions["G"] / (1 - base_proportions["C"]))  # type 2
    probG_to_C = base_proportions["G"] * (base_proportions["C"] / (1 - base_proportions["G"]))  # type 8
    probT_to_G = base_proportions["T"] * (base_proportions["G"] / (1 - base_proportions["T"]))  # type 6
    probA_to_C = base_proportions["A"] * (base_proportions["C"] / (1 - base_proportions["A"]))  # type 12
    class_probabilities["II"] = probC_to_G + probG_to_C + probT_to_G + probA_to_C

    probC_to_T = base_proportions["C"] * (base_proportions["T"] / (1 - base_proportions["C"]))  # type 3
    probG_to_A = base_proportions["G"] * (base_proportions["A"] / (1 - base_proportions["G"]))  # type 9
    probT_to_C = base_proportions["T"] * (base_proportions["C"] / (1 - base_proportions["T"]))  # type 5
    probA_to_G = base_proportions["A"] * (base_proportions["G"] / (1 - base_proportions["A"]))  # type 11
    class_probabilities["III"] = probC_to_T + probG_to_A + probT_to_C + probA_to_G

    print('\n' + "The probability of mutation classes I, II and III is: ", class_probabilities)
    return class_probabilities


def lr_class(type_counts: Dict[str, int], class_probabilities: Dict[str, float]):
    """Calculate the likelihood ratio of each mutation class - I,II,III.

    :param type_counts: dictionary where key is the mutation type and value is the total count of type
    :param class_probabilities: dictionary where key is class and value is the probability of the mutation class
    :return: lambda_mut_class, where key is mutation type and value is likelihood ratio of mutation class
    """
    # for each mutation class, tally the number of mutations
    total_I = total_II = total_III = 0
    for key in type_counts:
        if key in class_I_mutations:
            total_I += type_counts[key]
        elif key in class_II_mutations:
            total_II += type_counts[key]
        elif key in class_III_mutations:
            total_III += type_counts[key]

    # calculate likelihood ratio of each class
    sum_vtype = sum(type_counts.values())
    lambda_mut_class = {}
    for mut in mut_types:
        if mut in class_I_mutations:
            lambda_mut_class[mut] = total_I / (sum_vtype * class_probabilities["I"])
        elif mut in class_II_mutations:
            lambda_mut_class[mut] = total_II / (sum_vtype * class_probabilities["II"])
        elif mut in class_III_mutations:
            lambda_mut_class[mut] = total_III / (sum_vtype * class_probabilities["III"])
    print('\n' + "The likelihood of each mutation class is: ", lambda_mut_class)
    return lambda_mut_class


def count_type_per_pos(denovo_counts: Dict[str, int]):
    """Count the number of each mutation type at each mtDNA position - needed for make_mut_context_vector.
    The returned dictionary can be used for any region as the key includes the position.

    :param denovo_counts: a dictionary where the variant is the key, and the value is the total count of the variant
    :return: pos_by_type_counts, dictionary with count of each of the 12 mutation types at each mtDNA position
    """
    pos_by_type_counts = {}
    # use a list of all possible SNVs to parse quickly
    for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        variant = row["REF"] + row["POS"] + row["ALT"]
        mut_type = row["REF"] + ">" + row["ALT"]
        coord = row["POS"]
        # count the number of each 12 types of mutation at each position
        if variant in denovo_counts:
            allele_count = denovo_counts[variant]
            pos_by_type_counts[(coord, mut_type)] = allele_count
    return pos_by_type_counts


def make_mut_context_vector(pos_by_type_counts: Dict[Tuple[str, str], int],
                            flanking_range: List[int], reference_region: List[int], mut_group: str = None):
    """Generate the mutation sequence context vector Vseq t,p,n which counts the number of each flanking nucleotide
    around each mutation type for every flanking position, for both transitions and transversions.

    We measure sequence context separately for all four transitions, considering the strand of the reference nucleotide.
    For Tv, we calculate sequence context in a strand agnostic manner, for pyrimidine mutations at pyrimidine
    reference nucleotides C and T.

    :param pos_by_type_counts: dictionary with count of each of the 12 mutation types at each mtDNA position
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param reference_region: a list of the coordinates for the mtDNA region to calculate for
    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively, or the default of
    None which is to include both
    :return: vseq_vector, a dictionary with tuple key of mutation type, flanking nucleotide and flanking position,
    and value is the count for the tuple - for Tv the mutation keys are the pyrimidine mutations (C>A,C>G,T>A,T>G).
    """
    # first, generate dictionary to convert coordinate to reference nucleotide
    rcrs_pos2ref = rcrs_pos_to_ref()
    vseq_vector = {}
    for mut in mut_types:
        for coord in reference_region:
            if (str(coord), mut) in pos_by_type_counts:  # if the coordinate has a mutation
                for flanking_pos in flanking_range:
                    if mut in class_III_mutations:  # Ts
                        if mut_group == "Tv":
                            continue  # ie skip
                        # count number of A,C,G,T in the positions around each mutation type
                        # first, convert flanking position within [-10:10]/{0} to flanking coordinate
                        # then, convert the flanking coordinate to the flanking nucleotide in reference genome
                        # save flanking nucleotide counts around each mutation type, for each flanking position
                        flanking_coord = coord + flanking_pos
                        # this function handles positions near the artificial break and m.3107N spacer
                        flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                        flanking_nuc = rcrs_pos2ref[flanking_coord]  # reference nucleotide of flanking_position
                        vseq_vector = build_additive_dictionary(key=(mut, flanking_nuc, flanking_pos),
                                                                value=pos_by_type_counts[(str(coord), mut)],
                                                                dictionary=vseq_vector)
                    else:  # Tv
                        if mut_group == "Ts":
                            continue  # ie skip
                        # collapsed dictionary across strands for Tv, given lower counts
                        # generated as above, except annotate the mutations vs the pyrimidine (C or T)
                        if rcrs_pos2ref[str(coord)] == "C" or rcrs_pos2ref[str(coord)] == "T":
                            pyr_mut = mut
                            flanking_coord = coord + flanking_pos
                            flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                            flanking_nuc = rcrs_pos2ref[flanking_coord]
                            vseq_vector = build_additive_dictionary(key=(pyr_mut, flanking_nuc, flanking_pos),
                                                                    value=pos_by_type_counts[(str(coord), mut)],
                                                                    dictionary=vseq_vector)
                        elif rcrs_pos2ref[str(coord)] == "A" or rcrs_pos2ref[str(coord)] == "G":
                            pyr_mut = flip_mut(mut)  # reverse complement
                            flanking_coord = coord - flanking_pos  # reverse complement
                            flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                            flanking_nuc = flip_base(rcrs_pos2ref[flanking_coord])  # reverse complement
                            vseq_vector = build_additive_dictionary(key=(pyr_mut, flanking_nuc, flanking_pos),
                                                                    value=pos_by_type_counts[(str(coord), mut)],
                                                                    dictionary=vseq_vector)
    print('\n' + "The count of sequence context around mutations is: ", vseq_vector)
    return vseq_vector


def write_vseq_for_plotting(vseq_vector: Dict[Tuple[str, str, int], int], mut_type_counts: Dict[str, int],
                            flanking_range: List[int], region_name: str, mut_group: str = None):
    """Write the mutation sequence context vector Vseq t,p,n to file, for plotting.
    First convert from counts to frequency.

    :param vseq_vector: a dictionary with tuple key of mutation type, flanking nucleotide and flanking position,
    and value is the count for the tuple - for Tv the mutation key is the pyrimidine mutation (C>A,C>G,T>A,T>G)
    :param mut_type_counts: dictionary where key is the mutation type and the value is its total count in the region
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param region_name: the name of the reference region used, to include in file name
    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively, or the default of
    None which is to include both
    """
    # first do Ts, then Tv
    if mut_group is None or mut_group == "Ts":
        f = write_file_header(
            file='output_files/sequence_context_vectors/Vseq_mut_%s_%s.txt' % ("Ts", region_name),
            flanking_range=flanking_range, variable="mutation_type")
        for flanking_nuc in nucleotides:
            for mut in class_III_mutations:  # Ts
                f.write(flanking_nuc + '\t')
                for flanking_pos in flanking_range:
                    # to get as a frequency
                    vseq_freq = vseq_vector[(mut, flanking_nuc, flanking_pos)] / mut_type_counts[mut]
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(vseq_freq) + '\t' + mut + '\n')  # since end of row
                    else:
                        f.write(str(vseq_freq) + '\t')

    if mut_group is None or mut_group == "Tv":
        f = write_file_header(
            file='output_files/sequence_context_vectors/Vseq_mut_%s_%s.txt' % ("Tv", region_name),
            flanking_range=flanking_range, variable="mutation_type")
        for flanking_nuc in nucleotides:
            for mut in mut_types:
                if mut not in class_III_mutations:  # Tv
                    if mut[0] == "C" or mut[0] == "T":  # only iterate through pyrimidine Tv
                        f.write(flanking_nuc + '\t')
                        for flanking_pos in flanking_range:
                            # note vseq_vector only includes pyrimidine transversions (C>A, C>G, T>A, T>C)
                            # but mut_type_counts includes all 8 Tv types so sum the complements together (ie C>A+G>T)
                            vseq_freq = vseq_vector[(mut, flanking_nuc, flanking_pos)] / (
                                        mut_type_counts[mut] + mut_type_counts[flip_mut(mut)])
                            if flanking_pos == flanking_range[-1]:
                                f.write(str(vseq_freq) + '\t' + mut + '\n')  # since end of row
                            else:
                                f.write(str(vseq_freq) + '\t')


def make_ref_freq_vector(mut_group: str, flanking_range: List[int], reference_region: List[int]):
    """Generate the reference sequence context f(ref) n(t),p,n' which counts the number of each flanking nucleotide
    around each reference nucleotide for every flanking position.

    Do this separately for transitions - Ts - and transversions - Tv. For the Ts we measure sequence context separately
    for all four reference nucleotides A, C, G, T. For Tv, we calculate sequence context in a strand agnostic manner,
    for pyrimidine reference nucleotides C and T.

    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param reference_region: a list of the coordinates for the mtDNA region to count the occurrence of each base in
    :return: mut_group_fref_vector, a dictionary with tuple key of reference nucleotide, flanking nucleotide and
    flanking position, and value is the count for the tuple - for Tv the reference nucleotide is the pyrimidine (C or T)
    """
    # first, generate dictionary to convert coordinate to reference nucleotide
    rcrs_pos2ref = rcrs_pos_to_ref()
    mut_group_fref_vector = {}
    for coord in reference_region:
        if coord != 3107:  # where reference nucleotide is N
            ref_nuc = rcrs_pos2ref[str(coord)]
            for flanking_pos in flanking_range:
                if mut_group == "Ts":
                    # count number of A,C,G,T in the positions around each reference nucleotide A,C,G,T
                    # first, convert flanking position within [-10:10]/{0} to flanking coordinate
                    # then, convert the flanking coordinate to the flanking nucleotide in reference genome
                    # save flanking nucleotide counts around each reference nucleotide, for each flanking position
                    flanking_coord = coord + flanking_pos
                    # this function handles positions near the artificial break and m.3107N spacer
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    flanking_nuc = rcrs_pos2ref[flanking_coord]
                    mut_group_fref_vector = build_dictionary(key=(ref_nuc, flanking_nuc, flanking_pos),
                                                             dictionary=mut_group_fref_vector)
                elif mut_group == "Tv":
                    # collapsed dictionary across strands for Tv, given lower counts
                    # generated as above, except annotate the reference nucleotide vs the pyrimidine (C or T)
                    flanking_coord = pyr_ref_nuc = ''
                    if ref_nuc == "C" or ref_nuc == "T":
                        pyr_ref_nuc = ref_nuc
                        flanking_coord = coord + flanking_pos
                    elif ref_nuc == "A" or ref_nuc == "G":
                        pyr_ref_nuc = flip_base(ref_nuc)
                        flanking_coord = coord - flanking_pos  # reverse complement
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    flanking_nuc = rcrs_pos2ref[flanking_coord]
                    if ref_nuc == "A" or ref_nuc == "G":
                        flanking_nuc = flip_base(flanking_nuc)  # reverse complement
                    mut_group_fref_vector = build_dictionary(key=(pyr_ref_nuc, flanking_nuc, flanking_pos),
                                                             dictionary=mut_group_fref_vector)
    print('\n' + "The count of sequence context around reference nucleotides for mut group ", mut_group,
          " is: ", mut_group_fref_vector)
    return mut_group_fref_vector


def write_fref_for_plotting(mut_group: str, mut_group_fref_vector: Dict[Tuple[str, str, int], int],
                            base_proportions: Dict[str, float], flanking_range: List[int],
                            reference_region: List[int], region_name: str):
    """Write the reference sequence context f(ref) n(t),p,n' to file, for plotting.
    First convert from counts to frequency.

    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively
    :param mut_group_fref_vector: a dictionary with tuple key of reference nucleotide, flanking nucleotide and
    flanking position, and value is the count for the tuple - for Tv the reference nucleotide is the pyrimidine (C or T)
    :param base_proportions: a dictionary where key is the base and the value is the proportion of base in the region
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param reference_region: a list of the coordinates for the mtDNA region to count the occurrence of each base in
    :param region_name: the name of the reference region used, to include in file name
    """
    f = write_file_header(file='output_files/sequence_context_vectors/f_ref_%s_%s.txt' % (mut_group, region_name),
                          flanking_range=flanking_range, variable="reference_nucleotide")
    for flanking_nuc in nucleotides:
        if mut_group == "Ts":
            for ref_nuc in nucleotides:
                f.write(flanking_nuc + '\t')
                prob_ref_nuc = base_proportions[ref_nuc]
                for flanking_pos in flanking_range:
                    # to get as a frequency, before dividing by probability of reference nucleotide = A,C,G,T
                    f_ref = mut_group_fref_vector[(ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(f_ref / prob_ref_nuc) + '\t' + ref_nuc + '\n')  # since end of row
                    else:
                        f.write(str(f_ref / prob_ref_nuc) + '\t')
        elif mut_group == "Tv":
            for pyr_ref_nuc in ["C", "T"]:
                f.write(flanking_nuc + '\t')
                prob_ref_nuc = ''
                if pyr_ref_nuc == "C":
                    prob_ref_nuc = base_proportions["C"] + base_proportions["G"]
                elif pyr_ref_nuc == "T":
                    prob_ref_nuc = base_proportions["A"] + base_proportions["T"]
                for flanking_pos in flanking_range:
                    # to get as a frequency, before dividing by probability of pyrimidine reference nucleotide = C,T
                    f_ref = mut_group_fref_vector[(pyr_ref_nuc, flanking_nuc, flanking_pos)] / \
                            (len(reference_region) - 1)
                    if flanking_pos == flanking_range[-1]:
                        f.write(str(f_ref / prob_ref_nuc) + '\t' + pyr_ref_nuc + '\n')  # since end of row
                    else:
                        f.write(str(f_ref / prob_ref_nuc) + '\t')


def lr_seq_context(mut_type_counts: Dict[str, int], vseq_vector: Dict[Tuple[str, str, int], int],
                   Ts_fref_vector: Union[Dict[Tuple[str, str, int], int], None],
                   Tv_fref_vector: Union[Dict[Tuple[str, str, int], int], None],
                   base_proportions: Dict[str, float], flanking_range: List[int],
                   reference_region: List[int], mut_group: str = None):
    """Calculate the likelihood ratio of sequence context for all 12 mutation types.

    :param mut_type_counts: dictionary where key is the mutation type and the value is its total count in the region
    :param vseq_vector: dictionary with tuple key of mutation type, flanking nucleotide and flanking position,
    and value is the count for the tuple - for Tv the mutation is the pyrimidine mutation (C>A,C>G,T>A,T>G)
    :param Ts_fref_vector: a dictionary with tuple key of reference nucleotide, flanking nucleotide and
    flanking position, and value is the count for the tuple
    :param Tv_fref_vector: a dictionary with tuple key of reference nucleotide, flanking nucleotide and
    flanking position, and value is the count for the tuple - for Tv the reference nucleotide is the pyrimidine (C or T)
    :param base_proportions: a dictionary where key is the base and the value is the proportion of base in the region
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param reference_region: a list of the coordinates for the mtDNA region to count the occurrence of each base in
    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively, or the default of
    None which is to include both
    :return: lambda_seq_context, where key is a tuple of mutation type, flanking position and flanking nucleotide,
    and value is likelihood ratio of sequence context, and includes for all 12 mutation types
    """
    lambda_seq_context = {}
    for mut in mut_types:
        if mut in class_III_mutations:  # Ts
            if mut_group == "Tv":
                continue  # ie skip
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # numerator is sequence context for the base substitutions Vseq t,p,n
                    numerator = vseq_vector[(mut, flanking_nuc, flanking_pos)]
                    # denominator is Ts_fref_vector - the reference sequence context f(ref) n(t),p,n'
                    # divided by the number of positions to turn into frequency
                    # this is then divided by the probability of the reference nucleotide as a scaling factor
                    ref_nuc = mut[0]
                    f_ref = Ts_fref_vector[(ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    prob_ref = base_proportions[ref_nuc]
                    f_ref_scaled = f_ref / prob_ref
                    # and then multiplied by the Vtype count for the mutation type
                    denominator = mut_type_counts[mut] * f_ref_scaled
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = numerator / denominator
        elif mut in ["C>A", "C>G", "T>A", "T>G"]:  # pyrimidine Tv
            if mut_group == "Ts":
                continue  # ie skip
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # numerator is sequence context for the base substitutions Vseq t,p,n
                    # note vseq_vector includes all 4 Ts but only includes pyrimidine Tv (C>A, C>G, T>A, T>C)
                    numerator = vseq_vector[(mut, flanking_nuc, flanking_pos)]
                    # denominator is Tv_fref_vector - the reference sequence context f(ref) n(t),p,n'
                    # divided by the number of positions to turn into frequency
                    # this is then divided by the probability of the pyrimidine reference nucleotides as scaling factor
                    pyr_ref_nuc = mut[0]
                    f_ref = Tv_fref_vector[(pyr_ref_nuc, flanking_nuc, flanking_pos)] / (len(reference_region) - 1)
                    prob_ref = ''
                    # first need to collapse to the pyrimidine
                    if pyr_ref_nuc == "C":
                        prob_ref = base_proportions["C"] + base_proportions["G"]
                    elif pyr_ref_nuc == "T":
                        prob_ref = base_proportions["A"] + base_proportions["T"]
                    f_ref_scaled = f_ref / prob_ref
                    # and then multiplied by the Vtype count for the mutation type
                    # but mut_type_counts includes all 8 Tv types so sum the complements together (ie C>A+G>T)
                    denominator = (mut_type_counts[mut] + mut_type_counts[flip_mut(mut)]) * f_ref_scaled
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = numerator / denominator
    for mut in mut_types:  # have to wait until pyrimidine Tv are calculated
        if mut_group == "Ts":
            continue  # ie skip
        if mut in ["G>T", "G>C", "A>T", "A>C"]:  # purine Tv
            for flanking_pos in flanking_range:
                for flanking_nuc in nucleotides:
                    # these are the same values as for the reverse complement
                    lambda_seq_context[(mut, flanking_nuc, flanking_pos)] = \
                        lambda_seq_context[(flip_mut(mut), flip_base(flanking_nuc), (flanking_pos * -1))]
    print('\n' + "The likelihood of sequence context around mutation types is: ", lambda_seq_context)
    return lambda_seq_context


def write_lambda_seq_context_for_plotting(lambda_seq_context: Dict[Tuple[str, str, int], float],
                                          lambda_ref_nuc: Dict[str, float], lambda_mut_class: Dict[str, float],
                                          flanking_range: List[int], region_name: str,  mut_group: str = None):
    """Write the likelihood of sequence context to file, for plotting.
    Also include the likelihood of the reference nucleotide and likelihood of the mutation class.

    :param lambda_seq_context: dictionary where key is a tuple of mutation type, flanking nucleotide, flanking position,
    and value is likelihood ratio of sequence context
    :param lambda_ref_nuc: dictionary where key is reference nucleotide and value is likelihood ratio of mutation
    at the reference nucleotide
    :param lambda_mut_class: dictionary where key is mutation type and value is likelihood ratio of mutation class
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param region_name: the name of the reference region used, to include in file name
    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively, or the default of
    None which is to include both
    """
    f = write_file_header(
        file='output_files/sequence_context_vectors/lambda_seq_context_%s.txt' % region_name,
        flanking_range=flanking_range,
        variable="mutation_type" + '\t' + "lambda_ref_nucleotide" + '\t' + "lambda_mutation_class")
    for mut in mut_types:
        if mut_group == "Ts" and mut not in class_III_mutations:
            continue  # ie skip
        elif mut_group == "Tv" and mut in class_III_mutations:
            continue
        for flanking_nuc in nucleotides:
            f.write(str(flanking_nuc) + '\t')
            for flanking_pos in flanking_range:
                if flanking_pos == flanking_range[-1]:
                    f.write(str(lambda_seq_context[(mut, flanking_nuc, flanking_pos)]) + '\t' + str(mut) + '\t' +
                            str(lambda_ref_nuc[mut[0]]) + '\t' + str(lambda_mut_class[mut]) + '\n')
                else:
                    f.write(str(lambda_seq_context[(mut, flanking_nuc, flanking_pos)]) + '\t')


def inputs_for_composite_likelihood(denovo_list: str, reference_region: List[int], region_name: str,
                                    flanking_range: List[int], mut_group: str = None):
    """Function to make the dictionaries needed for the composite likelihood calculation.
    Configured to be able to run on just transitions, or transversions, or both which is the default.
    This is to handle the ori region separately.

    :param denovo_list: path to the final list of de novo mutations to use, in RefPosAlt format with their counts
    :param reference_region: a list of the coordinates for the mtDNA region to analyze
    :param region_name: the name of the reference region used, to include in file name
    :param flanking_range: a list of flanking positions to write data for, ranging up to [-10:10]/{0}
    :param mut_group: either Ts or Tv, to represent transitions or transversions respectively, or the default of
    None which is to include both
    :return: lambda_ref_nuc, lambda_mut_class, lambda_seq_context dictionaries for final function
    """
    print('Start with calculating likelihood of mutation at reference nucleotide and likelihood of mutation class')

    denovo_counts = make_denovo_counts(denovo_list=denovo_list)
    mut_type_counts = make_type_count_vector(denovo_counts=denovo_counts, reference_region=reference_region)
    base_proportions = probability_per_ref_nuc(reference_region=reference_region)
    lambda_ref_nuc = lr_ref_nuc(type_counts=mut_type_counts, base_proportions=base_proportions)
    class_probabilities = probability_per_class(base_proportions=base_proportions)
    lambda_mut_class = lr_class(type_counts=mut_type_counts, class_probabilities=class_probabilities)

    print('\n' + 'Now calculate likelihood of mutation sequence context')

    pos_by_type_counts = count_type_per_pos(denovo_counts=denovo_counts)
    vseq_vector = make_mut_context_vector(pos_by_type_counts=pos_by_type_counts, flanking_range=flanking_range,
                                          reference_region=reference_region, mut_group=mut_group)
    write_vseq_for_plotting(vseq_vector=vseq_vector, mut_type_counts=mut_type_counts, flanking_range=flanking_range,
                            region_name=region_name, mut_group=mut_group)

    # these functions have to have the mutation group specified, as handled separately
    Ts_fref_vector = Tv_fref_vector = {}
    if mut_group is None or mut_group == "Ts":
        Ts_fref_vector = make_ref_freq_vector(mut_group="Ts", flanking_range=flanking_range,
                                              reference_region=reference_region)
        write_fref_for_plotting(mut_group="Ts", mut_group_fref_vector=Ts_fref_vector,
                                base_proportions=base_proportions, flanking_range=flanking_range,
                                reference_region=reference_region, region_name=region_name)
    if mut_group is None or mut_group == "Tv":
        Tv_fref_vector = make_ref_freq_vector(mut_group="Tv", flanking_range=flanking_range,
                                              reference_region=reference_region)
        write_fref_for_plotting(mut_group="Tv", mut_group_fref_vector=Tv_fref_vector,
                                base_proportions=base_proportions, flanking_range=flanking_range,
                                reference_region=reference_region, region_name=region_name)

    lambda_seq_context = lr_seq_context(mut_type_counts=mut_type_counts, vseq_vector=vseq_vector,
                                        Ts_fref_vector=Ts_fref_vector, Tv_fref_vector=Tv_fref_vector,
                                        base_proportions=base_proportions, flanking_range=flanking_range,
                                        reference_region=reference_region, mut_group=mut_group)
    write_lambda_seq_context_for_plotting(lambda_seq_context=lambda_seq_context, lambda_ref_nuc=lambda_ref_nuc,
                                          lambda_mut_class=lambda_mut_class, flanking_range=flanking_range,
                                          region_name=region_name, mut_group=mut_group)

    return lambda_ref_nuc, lambda_mut_class, lambda_seq_context


def composite_likelihood(lambda_ref_nuc: Dict[str, float], lambda_mut_class: Dict[str, float],
                         lambda_seq_context: Dict[Tuple[str, str, int], float],
                         ori_lambda_ref_nuc: Dict[str, float], ori_lambda_mut_class: Dict[str, float],
                         ori_lambda_seq_context: Dict[Tuple[str, str, int], float],
                         context_size: int):
    """Calculate the composite likelihood of each mutation type at each position in the reference mtDNA sequence.

    :param lambda_ref_nuc: dictionary where key is reference nucleotide and value is likelihood ratio of mutation
    at the reference nucleotide
    :param lambda_mut_class: dictionary where key is mutation type and value is likelihood ratio of mutation class
    :param lambda_seq_context: where key is a tuple of mutation type, flanking nucleotide, and flanking position,
    and value is likelihood ratio of sequence context, and includes for all 12 mutation types
    :param ori_lambda_ref_nuc: dictionary where key is reference nucleotide and value is likelihood ratio of mutation
    at the reference nucleotide, in ori region only
    :param ori_lambda_mut_class: dictionary where key is mutation type and value is likelihood ratio of mutation class,
    in ori region only
    :param ori_lambda_seq_context: where key is a tuple of mutation type, flanking nucleotide, and flanking position,
    and value is likelihood ratio of sequence context, in ori region only
    :param context_size: how many nucleotides to include for sequence context, 3 is default (trinucleotide)
    """
    f = open('output_files/mutation_likelihoods/mito_mutation_likelihoods.txt', "w")
    f.write("POS	REF	ALT	Likelihood" + '\n')

    # how many nucleotides either side to include in model, +1 since range end not included
    pos_range = list(range(1, ((context_size - 1) // 2) + 1))

    # first, generate dictionary to convert coordinate to reference nucleotide
    rcrs_pos2ref = rcrs_pos_to_ref()

    # iterate through all possible mutations in the mitochondrial genome
    for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        coord = int(row["POS"])
        ref = row["REF"]
        alt = row["ALT"]
        mut = ref + ">" + alt

        if coord != 3107:  # this is ref="N"
            # product_mut_type is likelihood of mutation at reference nucleotide
            # multiplied by likelihood of mutation class - ori region handled separately
            if coord in ori_region:
                product_mut_type = ori_lambda_ref_nuc[ref] * ori_lambda_mut_class[mut]
            else:
                product_mut_type = lambda_ref_nuc[ref] * lambda_mut_class[mut]

            # reset product_seq_context to 0 for each mutation type at each reference coordinate
            product_seq_context = 0

            for number in pos_range:  # this list will range from 1-10, depending on size of sequence context
                # convert the flanking positions to their coordinates
                flanking_coord_list = [coord - number, coord + number]
                for flanking_coord in flanking_coord_list:
                    # to designate which is the minus and plus flanking position around the reference coordinate
                    flanking_pos = ''
                    if flanking_coord < coord:
                        flanking_pos = -1 * number
                    elif flanking_coord > coord:
                        flanking_pos = number
                    # this function handles positions near the artificial break and m.3107N spacer
                    flanking_coord = handle_mito_genome(flanking_coord=flanking_coord, coord=coord)
                    # reference nucleotide of flanking_position
                    flanking_nuc = rcrs_pos2ref[flanking_coord]

                    # calculate product of the likelihood for sequence context for each flanking position
                    # Ts in ori region handled separately
                    if (coord in ori_region) and (mut in class_III_mutations):  # if Ts in OriB-OriH
                        if product_seq_context == 0:
                            product_seq_context = ori_lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                        else:
                            product_seq_context = product_seq_context * \
                                                  ori_lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                    else:
                        if product_seq_context == 0:
                            product_seq_context = lambda_seq_context[(mut, flanking_nuc, flanking_pos)]
                        else:
                            product_seq_context = product_seq_context * \
                                                  lambda_seq_context[(mut, flanking_nuc, flanking_pos)]

            # final composite likelihood - likelihood mutation type multiplied by likelihood sequence context
            product_pos = product_mut_type * product_seq_context
            # print to file:
            f.write(str(coord) + '\t' + str(ref) + '\t' + str(alt) + '\t' + str(product_pos) + '\n')


if __name__ == "__main__":
    print(datetime.datetime.now(), '\n' + "Starting to build composite likelihood model for mtDNA!" + '\n')

    parser = argparse.ArgumentParser()
    parser.add_argument("-context_size", type=int,
                        help="how many nucleotides to include when calculating sequence context, ie 3 is trinucleotide")
    parser.add_argument("-denovo_list", type=str,
                        help="path to the de novo mutation list, in RefPosAlt format with their counts")
    args = parser.parse_args()

    # set defaults
    if args.context_size is None:
        args.context_size = 3
    if args.denovo_list is None:
        args.denovo_list = 'output_files/denovo/final_denovo.txt'

    for path in ['output_files/sequence_context_vectors/', 'output_files/mutation_likelihoods/']:
        if not os.path.exists(path):
            os.makedirs(path)
    print("Creating required directories")

    print('\n' + 'This script will calculate the likelihood of mutation in the mitochondrial genome')
    print('\n' + 'This will be done using window size of ' + str(args.context_size) + ' for sequence context')
    # double bracket division to return as int
    flanking_range = [i for i in
                      list(range(-((args.context_size - 1) // 2), (((args.context_size - 1) // 2) + 1)))
                      if i != 0]
    print('\n' + 'Mutation likelihoods in the OriB-OriH region from ' + str(start_ori) + '-' + str(end_ori) +
          ' will be calculated separately')

    print('\n' + 'First, start with the reference region excluding the OriB-OriH' + '\n')

    (lambda_ref_nuc, lambda_mut_class, lambda_seq_context) = \
        inputs_for_composite_likelihood(denovo_list=args.denovo_list, reference_region=reference_except_ori,
                                        region_name="", flanking_range=flanking_range)

    print('\n' + 'Next, calculate for the OriB-OriH region' + '\n')

    # only compute this for transitions - Ts
    (ori_lambda_ref_nuc, ori_lambda_mut_class, ori_lambda_seq_context) = \
        inputs_for_composite_likelihood(denovo_list=args.denovo_list, reference_region=ori_region,
                                        region_name="OriB-OriH", flanking_range=flanking_range, mut_group="Ts")

    print('\n' + 'Now, calculate position mutability across mitochondrial genome and write to file' + '\n')

    composite_likelihood(lambda_ref_nuc=lambda_ref_nuc, lambda_mut_class=lambda_mut_class,
                         lambda_seq_context=lambda_seq_context, ori_lambda_ref_nuc=ori_lambda_ref_nuc,
                         ori_lambda_mut_class=ori_lambda_mut_class, ori_lambda_seq_context=ori_lambda_seq_context,
                         context_size=args.context_size)

    print(datetime.datetime.now(), '\n' + "Finished building composite likelihood model for mtDNA!")
