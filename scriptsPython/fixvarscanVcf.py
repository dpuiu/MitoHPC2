#!/usr/bin/env python3

import sys
import pysam

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():
    """
    Filter VCF records from stdin using pysam.

    Filters applied:
    - Removes variants at specific positions (3105, 3106)
    - Removes records where REF contains 'N'

    Input:
        VCF streamed via stdin (pysam.VariantFile)

    Output:
        Filtered VCF to stdout with original header preserved

    Example:
        cat I.vcf | python script.py > O.vcf
    """

    try:
        vcf_in = pysam.VariantFile(sys.stdin)
        vcf_out = pysam.VariantFile(sys.stdout, "w", header=vcf_in.header)

        for rec in vcf_in:
            # Skip specific positions (VCF POS is 1-based)
            if rec.pos in (3105, 3106):
                continue

            # Skip records where REF contains 'N'
            if rec.ref and "N" in rec.ref:
                continue

            vcf_out.write(rec)

        vcf_out.close()

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
