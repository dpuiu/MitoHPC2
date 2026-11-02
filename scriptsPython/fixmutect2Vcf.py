#!/usr/bin/env python3

import argparse
import signal
from typing import Dict, Any

import pysam

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


KEEP: tuple[str, ...] = ("GT", "DP", "AD", "AF")


def build_sample(sample: pysam.libcbcf.VariantRecordSample) -> Dict[str, Any]:
    """
    Extract only required Mutect2 fields.
    """

    return {
        "GT": sample.get("GT"),
        "DP": sample.get("DP"),
        "AD": sample.get("AD"),
        "AF": sample.get("AF"),
    }


def process_vcf(
    vcf_in: pysam.VariantFile,
    vcf_out: pysam.VariantFile,
) -> None:
    """
    Stream VCF and keep only GT:DP:AD:AF FORMAT fields.
    """

    for rec in vcf_in:

        if rec.alts is None:
            continue

        if len(rec.samples) == 0:
            continue

        sample = rec.samples[0]

        fixed: Dict[str, Any] = build_sample(sample)

        # REMOVE ALL FORMAT FIELDS FIRST (critical fix)
        for k in list(sample.keys()):
            del sample[k]

        # REINSERT ONLY REQUIRED ONES
        for k in KEEP:
            sample[k] = fixed[k]


        vcf_out.write(rec)


def main() -> None:
    """
    CLI entry point.
    """

    parser = argparse.ArgumentParser(
        description="Mutect2 VCF cleaner (GT:DP:AD:AF only)"
    )

    parser.add_argument("-i", "--input", default="-")
    parser.add_argument("-o", "--output", default="-")

    args = parser.parse_args()

    vcf_in = pysam.VariantFile(args.input, "r")
    vcf_out = pysam.VariantFile(args.output, "w", header=vcf_in.header)

    process_vcf(vcf_in, vcf_out)

    vcf_out.close()


if __name__ == "__main__":
    main()
