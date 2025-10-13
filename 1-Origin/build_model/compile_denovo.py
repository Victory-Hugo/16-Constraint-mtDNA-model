import csv
import datetime
from typing import Dict, TextIO, Union
import os


# helper functions
def build_sample_dict(sample: str, dict: Dict[str, int]):
    """Generate dictionary with the de novo count of each sample.

    :param sample: sample name
    :param dict: dictionary name
    :return: a dictionary, where the key is the sample, and the value is the sample's de novo count
    """
    if sample not in dict:
        dict[sample] = 1
    else:
        dict[sample] += 1
    return dict


def write_denovo(file: TextIO, sample: str, variant: str, dict: Union[Dict[str, int], int]):
    """Writes out de novo variants alongside sample it was found in, and the de novo count of that sample.

    :param file: file that the de novo variants will be written to
    :param sample: sample name
    :param variant: de novo variant, in RefPosAlt format
    :param dict: a dictionary, where the key is the sample, and the value is the sample's de novo count - the count can
    also be directly provided
    """
    if type(dict) is int:
        sample_count = dict
    else:
        sample_count = dict[sample]
    file.write(variant + '\t' + sample + '\t' + str(sample_count) + '\n')


def rcrs_pos_to_ref():
    """Generate dictionary linking each position to its reference nucleotide in the rCRS.

    :return: dictionary where the key is the position in rCRS, and the value is its reference nucleotide
    """
    dictionary = {}
    for row in csv.DictReader(open('required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf'),
                              delimiter='\t'):
        dictionary[row["POS"]] = row["REF"]
    return dictionary


# function for compiling de novo
def extract_denovo(file: TextIO):
    """Function to parse the de novo lists obtained from the literature, and own datasets.
    Each file parsed uniquely, given they are all in different formats.
    Three categories - germline, somatic tissue, and somatic cancer de novo variants.

    :param file: the txt file to write the de novo variants, their sample, and the de novo count of that sample
    """
    # first, build the sample_counts dictionary with the count of de novos for each sample
    # then, write the de novo variant and the sample details to file
    sample_counts = {}

    # GERMLINE

    # dataset 1 - from Wei et al 2019 Science

    for row in csv.reader(open('required_files/input_denovo/germline/PMID31123110_Data_S1.txt'), delimiter='\t'):
        if row[6] == "de novo":
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID31123110" + "-germline"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID31123110_Data_S1.txt'), delimiter='\t'):
        if row[6] == "de novo":
            write_denovo(file=file, sample=(row[0] + "-PMID31123110" + "-germline"), variant=row[2], dict=sample_counts)

    # dataset 2 - from Rebolledo-Jaramillo et al 2014 PNAS

    for row in csv.reader(open('required_files/input_denovo/germline/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "child":  # category child = germline de novo variants, see manuscript
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID25313049" + "-germline"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "child":
            write_denovo(file=file, sample=(row[0] + "-PMID25313049" + "-germline"), variant=(row[2] + row[1] + row[3]),
                         dict=sample_counts)

    # dataset 3 - from Zaidi et al 2019 PNAS

    catch_list = []  # use this to catch the de novo variants listed twice for one individual (one line per tissue type)
    for row in csv.reader(open('required_files/input_denovo/germline/PMID31757848_TableS3.txt'), delimiter='\t'):
        if row[6] == "1":  # heuristic method of identification, ie from sequencing, see manuscript
            sample = row[2] + "-PMID31757848" + "-germline"
            variant = row[9] + row[4] + row[10]
            if (sample, variant) not in catch_list:
                catch_list.append((sample, variant))
                sample_counts = build_sample_dict(sample=sample, dict=sample_counts)

    catch_list = []
    for row in csv.reader(open('required_files/input_denovo/germline/PMID31757848_TableS3.txt'), delimiter='\t'):
        if row[6] == "1":
            sample = row[2] + "-PMID31757848" + "-germline"
            variant = row[9] + row[4] + row[10]
            if (sample, variant) not in catch_list:
                catch_list.append((sample, variant))
                write_denovo(file=file, sample=sample, variant=variant, dict=sample_counts)

    # dataset 4 - from Li et al 2016 Gen Res

    for row in csv.reader(open('required_files/input_denovo/germline/PMID26916109_Dataset_S1_Heteroplasmy_list.txt'),
                          delimiter='\t'):
        if row[8] == "Yes":  # de novo
            sample_counts = build_sample_dict(sample=(row[0] + row[1] + "-PMID26916109" + "-germline"),
                                              dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/germline/PMID26916109_Dataset_S1_Heteroplasmy_list.txt'),
                          delimiter='\t'):
        if row[8] == "Yes":  # de novo
            write_denovo(file=file, sample=(row[0] + row[1] + "-PMID26916109" + "-germline"),
                         variant=(row[5] + row[3] + row[6]), dict=sample_counts)

    # dataset 5 - from my own SPARK analyses
    # note, this cannot be released due to data restrictions
    if os.path.exists('required_files/input_denovo/germline/SPARK_mtDNA_de_novo.txt'):
        for row in csv.DictReader(open('required_files/input_denovo/germline/SPARK_mtDNA_de_novo.txt'), delimiter='\t'):
            s_counts = int(row["number"])  # taken directly rather than within dictionary
            write_denovo(
                file=file, sample=("sample_from" + "-SPARK" + "-germline"), variant=row["de_novo"], dict=s_counts)
        
        
    # SOMATIC TISSUE

    # dataset 2 from germline (from Rebolledo-Jaramillo et al 2014 PNAS)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "somatic-gain":  # category somatic-gain = somatic de novo variants, see manuscript
            sample_counts = build_sample_dict(sample=(row[0] + "-PMID25313049" + "-somatic_tissue"), dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID25313049_TableS3.txt'), delimiter='\t'):
        if row[10] == "somatic-gain":
            write_denovo(file=file, sample=(row[0] + "-PMID25313049" + "-somatic_tissue"),
                         variant=(row[2] + row[1] + row[3]), dict=sample_counts)

    # dataset 3 from germline (from Zaidi et al 2019 PNAS)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID31757848_TableS4.txt'), delimiter='\t'):
        sample_counts = build_sample_dict(sample=(row[0] + row[4] + "-PMID31757848" + "-somatic_tissue"),
                                          dict=sample_counts)

    for row in csv.reader(open('required_files/input_denovo/somatic_tissue/PMID31757848_TableS4.txt'), delimiter='\t'):
        if not (row[0].startswith('Table S4')) and not (row[0].startswith('SampleID')):  # to skip header lines
            write_denovo(file=file, sample=(row[0] + row[4] + "-PMID31757848" + "-somatic_tissue"),
                         variant=(row[7] + row[5] + row[8]), dict=sample_counts)

    # dataset 6 - gtex, from Ludwig et al 2019 Cell

    for row in csv.reader(
            open('required_files/input_denovo/somatic_tissue/PMID30827679_NIHMS1518665-supplement-10.txt'),
            delimiter='\t'):
        if not (row[0].startswith('Table')) and not (row[0].startswith('Mutation')):
            sample_counts = build_sample_dict(sample=(row[1] + "-PMID30827679" + "-somatic_tissue"), dict=sample_counts)

    # need rCRS lookup for the conversion of de novo variants from the GTEx dataset from Yoruba reference to rCRS
    rcrs_pos2ref = rcrs_pos_to_ref()
    # positions with differences between Yoruba Sequence (L3e2b1a1) to rCRS (H2a2a1)
    # per https://haplogrep.i-med.ac.at/2014/09/08/rcrs-vs-rsrs-vs-hg19/
    yoruba_rCRS_diff = [73, 150, 195, 263, 309, 315, 408, 750, 1438, 2352, 2483, 2706, 3107, 4769, 5580, 7028, 8701,
                        8860, 9377, 9540, 10398, 10819, 10873, 11017, 11719, 11722, 12705, 12850, 14212, 14580, 14766,
                        14905, 15301, 15326, 15932, 16172, 16183, 16189, 16193, 16223, 16320, 16519]

    for row in csv.reader(
            open('required_files/input_denovo/somatic_tissue/PMID30827679_NIHMS1518665-supplement-10.txt'),
            delimiter='\t'):
        if not (row[0].startswith('Table')) and not (row[0].startswith('Mutation')):
            pos = int(row[0].split('_')[0])
            alt = row[0].split('_')[1]

            # appears to be the Yoruban reference sequence, so convert to rCRS
            # per https://www.mitomap.org/foswiki/bin/view/MITOMAP/YorubanConversion
            if 311 <= pos <= 316:  # before 309 is the same
                pos = pos - 1
            elif 318 <= pos <= 3108:
                pos = pos - 2
            elif 3109 <= pos <= 16190:
                pos = pos - 1
            elif 16192 <= pos <= 16571:
                pos = pos - 2
            elif pos == 310 or pos == 317 or pos == 16191:
                continue
            if pos in yoruba_rCRS_diff:
                continue  # want to skip these

            ref = rcrs_pos2ref[str(pos)]
            write_denovo(file=file, sample=(row[1] + "-PMID30827679" + "-somatic_tissue"),
                         variant=(ref + str(pos) + alt), dict=sample_counts)

    # SOMATIC CANCER

    # dataset 7 - from Yuan et al 2020 Nature Genetics

    for row in csv.DictReader(open('required_files/input_denovo/somatic_cancer/TCMA-MutationSNV.tsv'), delimiter='\t'):
        sample_counts = build_sample_dict(sample=(row["sample_id"] + '-PMID32024997' + '-somatic_cancer'),
                                          dict=sample_counts)

    for row in csv.DictReader(open('required_files/input_denovo/somatic_cancer/TCMA-MutationSNV.tsv'), delimiter='\t'):
        sample = row["sample_id"] + '-PMID32024997' + '-somatic_cancer'
        # remove 7 samples authors note to be hypermutated, defined as >13 somatic variants
        if sample_counts[sample] < 14:
            write_denovo(file=file, sample=sample, variant=(row["ref"] + row["position"] + row["var"]),
                         dict=sample_counts)


if __name__ == "__main__":
    print(datetime.datetime.now(), "Starting to compile all de novo!")

    if not os.path.exists('output_files/denovo'):
        os.makedirs('output_files/denovo')
    print(datetime.datetime.now(), "Creating required directories")

    # create file to write all de novo variants and their sample details
    f = open('output_files/denovo/all_denovo.txt', "w")
    f.write("denovo	sample	sample_denovo_count" + '\n')

    extract_denovo(file=f)

    print(datetime.datetime.now(), "Script finished!")
