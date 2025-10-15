import argparse
import collections
import csv
import gzip
import os
from typing import Dict, List, Tuple


VariantCarrierStats = collections.namedtuple(
    "VariantCarrierStats", ["carrier_count", "callable_samples", "allele_frequency"]
)


def rcrs_pos_to_ref() -> Dict[str, str]:
    """Build a lookup table mapping rCRS positions to reference bases."""

    lookup: Dict[str, str] = {}
    with open('0-required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            lookup[row["POS"]] = row["REF"]
    return lookup


def ensure_output_dirs():
    """Create the output directories expected by downstream scripts."""
    for path in ["output", "output/denovo"]:
        if not os.path.exists(path):
            os.makedirs(path)


def parse_vcf(
    vcf_path: str,
    max_af: float,
    min_carriers: int,
) -> Tuple[Dict[str, VariantCarrierStats], Dict[str, List[str]]]:
    """Parse a multi-sample mtDNA VCF and collect carrier information per variant and sample.

    Parameters
    ----------
    vcf_path : str
        Path to the bgzip-compressed VCF containing mtDNA variation.
    max_af : float
        Maximum allele frequency threshold (AC / AN) for a variant to be retained.
    min_carriers : int
        Minimum number of carriers required for a variant to be retained.

    Returns
    -------
    Tuple containing:
        1. Dictionary keyed by variant string (RefPosAlt) with carrier statistics.
        2. Dictionary keyed by sample with a list of retained variant strings.
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    variant_stats: Dict[str, VariantCarrierStats] = {}
    sample_variants: Dict[str, List[str]] = collections.defaultdict(list)

    with gzip.open(vcf_path, "rt") as handle:
        samples: List[str] = []
        gt_index = None

        for raw_line in handle:
            line = raw_line.rstrip()
            if not line or line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.split("\t")
                samples = header[9:]
                for sample in samples:
                    sample_variants[sample] = []
                continue

            fields = line.split("\t")
            chrom, pos, _id, ref, alt, _qual, flt, info, fmt = fields[:9]
            genotypes = fields[9:]

            if chrom not in {"chrM", "MT", "M"}:
                continue
            if flt not in {"PASS", "."}:
                continue

            if "," in alt:
                # Skip multi-allelic sites for now.
                continue

            if len(ref) != 1 or len(alt) != 1:
                # Restrict to SNVs for the current carrier-count model.
                continue

            info_dict = {item.split("=")[0]: item.split("=")[1] for item in info.split(";") if "=" in item}
            ac = int(info_dict.get("AC", 0))
            an = int(info_dict.get("AN", 0))

            if an == 0 or ac == 0:
                continue

            af = ac / an
            if af > max_af or ac < min_carriers:
                continue

            format_keys = fmt.split(":")
            if gt_index is None:
                try:
                    gt_index = format_keys.index("GT")
                except ValueError as exc:
                    raise ValueError("FORMAT column does not define GT field") from exc

            variant_key = (ref, pos, alt)
            variant_id = f"{ref}{pos}{alt}"

            callable_samples = 0
            carrier_count = 0

            for sample, genotype_str in zip(samples, genotypes):
                genotype_fields = genotype_str.split(":")
                gt = genotype_fields[gt_index]

                if gt in {".", "./.", ".|."}:
                    continue

                callable_samples += 1

                alleles = gt.replace("|", "/").split("/")
                if any(allele not in {"0", "."} for allele in alleles):
                    carrier_count += 1
                    sample_variants[sample].append(variant_id)

            if carrier_count == 0 or callable_samples == 0:
                continue

            allele_frequency = carrier_count / callable_samples
            variant_stats[variant_key] = VariantCarrierStats(
                carrier_count=carrier_count,
                callable_samples=callable_samples,
                allele_frequency=allele_frequency,
            )

    return variant_stats, sample_variants


def write_all_denovo(sample_variants: Dict[str, List[str]], output_path: str) -> None:
    """Write per-sample variant records in the legacy all_denovo format."""
    with open(output_path, "w") as handle:
        handle.write("denovo\tsample\tsample_denovo_count\n")
        for sample, variants in sample_variants.items():
            if not variants:
                continue
            sample_variant_count = len(variants)
            for variant in variants:
                handle.write(f"{variant}\t{sample}\t{sample_variant_count}\n")


def write_final_denovo(variant_stats: Dict[str, VariantCarrierStats], output_path: str) -> None:
    """Write the aggregated carrier counts per variant."""
    with open(output_path, "w") as handle:
        handle.write("denovo\tcount\n")
        for (ref, pos, alt), stats in sorted(variant_stats.items(), key=lambda item: int(item[0][1])):
            variant = f"{ref}{pos}{alt}"
            handle.write(f"{variant}\t{stats.carrier_count}\n")


def write_variant_stats(variant_stats: Dict[str, VariantCarrierStats], output_path: str) -> None:
    """Persist variant statistics for downstream annotation."""
    with open(output_path, "w") as handle:
        handle.write("POS\tREF\tALT\tcarrier_count\tcallable_samples\tallele_frequency\n")
        for (ref, pos, alt), stats in sorted(variant_stats.items(), key=lambda item: int(item[0][1])):
            handle.write(
                f"{pos}\t{ref}\t{alt}\t{stats.carrier_count}\t"
                f"{stats.callable_samples}\t{stats.allele_frequency}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description="Convert a multi-sample mtDNA VCF into carrier-count inputs for the constraint pipeline."
    )
    parser.add_argument(
        "--vcf",
        dest="vcf_path",
        default="example.vcf.gz",
        help="Path to the bgzip-compressed mtDNA VCF (default: example.vcf.gz in the project root).",
    )
    parser.add_argument(
        "--max-af",
        dest="max_af",
        type=float,
        default=0.01,
        help="Maximum allele frequency (AC/AN) for variants to be retained (default: 0.01).",
    )
    parser.add_argument(
        "--min-carriers",
        dest="min_carriers",
        type=int,
        default=1,
        help="Minimum number of carriers required to retain a variant (default: 1).",
    )
    args = parser.parse_args()

    ensure_output_dirs()

    print("Parsing VCF and collecting carrier information...")
    variant_stats, sample_variants = parse_vcf(
        vcf_path=args.vcf_path,
        max_af=args.max_af,
        min_carriers=args.min_carriers,
    )

    print(f"Retained variants: {len(variant_stats)}")

    write_all_denovo(sample_variants, "output/denovo/all_denovo.txt")
    write_final_denovo(variant_stats, "output/denovo/final_denovo.txt")
    write_variant_stats(variant_stats, "output/denovo/variant_carrier_stats.txt")

    print("Carrier-based inputs have been written to output/denovo/")


if __name__ == "__main__":
    main()
