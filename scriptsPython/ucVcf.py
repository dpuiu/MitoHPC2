#!/usr/bin/env python3

import sys
import pysam

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():
    try:
        vcf_in = pysam.VariantFile(sys.stdin)
        vcf_out = pysam.VariantFile(sys.stdout, "w", header=vcf_in.header)

        for rec in vcf_in:
            rec.ref = rec.ref.upper()
            if rec.alts:
                rec.alts = tuple(a.upper() for a in rec.alts)
            vcf_out.write(rec)

        vcf_out.close()

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
