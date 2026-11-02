#!/usr/bin/env python3

import sys
import pysam
import re

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():
    """
    Reformat and filter a VCF stream using pysam.

    Functionality:
    - Removes variants where REF contains 'N'
    - Validates FORMAT field order
    - Reorders FORMAT fields to: GT:DP:AD:AF
    - Rewrites sample fields accordingly
    - Writes updated VCF to stdout with original header

    Input:
        VCF from stdin (pysam.VariantFile)

    Output:
        Modified VCF stream to stdout

    Expected FORMAT order before transformation:
        GT:GQ:DP:AF:AD

    Output FORMAT order:
        GT:DP:AD:AF

    Raises:
        ValueError if FORMAT structure is not as expected

    Example:
        cat I.vcf | python script.py > O.vcf
    """
    try:
        vcf_in = pysam.VariantFile(sys.stdin)
        vcf_out = pysam.VariantFile(sys.stdout, "w", header=vcf_in.header)

        for rec in vcf_in:

            # skip REF containing N
            if rec.ref and "N" in rec.ref:
                continue

            # ensure expected FORMAT order exists
            fmt_keys = list(rec.format.keys()) if rec.format else []
            if fmt_keys[:5] != ["GT", "GQ", "DP", "AF", "AD"]:
                # equivalent to Perl "die"
                raise ValueError(f"Unexpected FORMAT header in record: {rec}")

            sample = list(rec.samples.values())[0]

            gt = sample.get("GT")
            dp = sample.get("DP")
            af = sample.get("AF")
            ad = sample.get("AD")

            # rebuild FORMAT like Perl script
            rec.format = ["GT", "DP", "AD", "AF"]

            # update sample fields in new order
            sample.clear()
            sample["GT"] = gt
            sample["DP"] = dp
            sample["AD"] = ad
            sample["AF"] = af

            vcf_out.write(rec)

        vcf_out.close()

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
