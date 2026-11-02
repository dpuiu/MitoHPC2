#!/usr/bin/env python3

import argparse
import signal
from typing import Any, Dict

import pysam

signal.signal(signal.SIGPIPE, signal.SIG_DFL)


# -------------------------
# Helpers
# -------------------------
def extract_af(sample: Any) -> float:
    """
    Extract allele frequency (AF) from a pysam sample.

    Args:
        sample: pysam sample object (untyped)

    Returns:
        AF value (float), default = 1.0 if missing
    """
    try:
        af = sample.get("AF")
        return float(af) if af is not None else 1.0
    except Exception:
        return 1.0


# -------------------------
# Core logic
# -------------------------
def process(vcf_in: pysam.VariantFile) -> None:
    """
    Filter VCF records by allele frequency dominance.

    Logic:
        - track max AF per position
        - track sum AF per position
        - keep if: max_AF > (1 - sum_AF)
        - enforce non-overlapping output ordering
    """

    max_af: Dict[int, float] = {}
    sum_af: Dict[int, float] = {}
    best_rec: Dict[int, Any] = {}

    # pass 1: collect statistics
    for rec in vcf_in:

        if rec.alts is None:
            continue

        pos: int = rec.pos
        sample = rec.samples[0]
        af = extract_af(sample)

        if pos not in max_af or af > max_af[pos]:
            max_af[pos] = af
            best_rec[pos] = rec

        sum_af[pos] = sum_af.get(pos, 0.0) + af

    # filter step
    selected: Dict[int, Any] = {}

    for pos in best_rec:
        if max_af[pos] > (1.0 - sum_af[pos]):
            selected[pos] = best_rec[pos]

    # output with overlap control
    pos_max = 0

    for pos in sorted(selected):

        rec = selected[pos]

        ref_len = len(rec.ref)
        alt_len = len(rec.alts[0]) if rec.alts else 1

        diff = abs(ref_len - alt_len)

        if pos >= pos_max:
            print(
                rec.chrom,
                rec.pos + 1,
                rec.id,
                rec.ref,
                rec.alts[0],
                rec.qual if rec.qual is not None else ".",
                sep="\t",
            )

        pos_max = max(pos_max, pos + diff)


# -------------------------
# CLI
# -------------------------
def main() -> None:
    """
    Entry point:
    - parses arguments
    - streams VCF
    - applies AF filter
    """

    parser = argparse.ArgumentParser(
        description="AF-based VCF filter (Mutect2-style)"
    )

    parser.add_argument("-i", "--input", default="-")

    args = parser.parse_args()

    vcf_in = pysam.VariantFile(args.input, "r")

    process(vcf_in)


if __name__ == "__main__":
    main()
