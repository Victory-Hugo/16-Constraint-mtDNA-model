import argparse
from compare_denovo import filter_denovo
import datetime


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-germline_max", type=int,
                        help="Maximum sample de novo count for inclusion, for germline de novo")
    parser.add_argument("-som_tissue_max", type=int,
                        help="Maximum sample de novo count for inclusion, for somatic tissue de novo")
    parser.add_argument("-som_cancer_max", type=int,
                        help="Maximum sample de novo count for inclusion, for somatic cancer de novo")
    args = parser.parse_args()

    # set default, note germline_max and som_tissue_max are already None and don't want to apply a threshold for them
    if args.som_cancer_max is None:
        args.som_cancer_max = 1

    print(datetime.datetime.now(), "Starting to filter de novo!")
    print('\n' + "This will produce a final list of de novo that are used to calculate mutational likelihood scores")

    f = open('output_files/denovo/final_denovo.txt', "w")
    f.write("denovo	count" + '\n')

    denovo_counts = {}
    denovo_counts = filter_denovo(germline_max=args.germline_max, som_tissue_max=args.som_tissue_max,
                                  som_cancer_max=args.som_cancer_max, denovo_counts=denovo_counts)

    for key in denovo_counts:
        f.write(key + '\t' + str(denovo_counts[key]) + '\n')

    print(datetime.datetime.now(), "Script finished!")
