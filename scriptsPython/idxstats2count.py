#!/usr/bin/env python3

import sys
import argparse
import os

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():
    """
    Parse a samtools idxstats file and compute:

    - total reads (mapped + unmapped across all contigs)
    - total mapped reads
    - mapped reads for a specific chromosome (default: chrM)

    Input
    -----
    Reads from stdin in samtools idxstats format:
        chrom \t length \t mapped \t unmapped

    Output
    ------
    Tab-delimited summary:
        Run    all_reads    mapped_reads    MT_reads

    Example
    -------
    cat I.idxstats | python {script_name} -sample S1 --chrM chrM
    """.format(script_name=os.path.basename(__file__))

    parser = argparse.ArgumentParser(
        description="Count total, mapped, and chrM reads from samtools idxstats"
    )

    parser.add_argument("--sample", default=".")
    parser.add_argument("--chrM", default="chrM")

    args = parser.parse_args()

    total_reads = 0
    mapped_reads = 0
    chrM_reads = 0

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        fields = line.split("\t")
        if len(fields) < 4:
            continue

        chrom = fields[0]
        mapped = int(fields[2])
        unmapped = int(fields[3])

        total_reads += mapped + unmapped
        mapped_reads += mapped

        if chrom == args.chrM:
            chrM_reads = mapped

    # output (same format as Perl)
    print("\t".join(["Run", "all_reads", "mapped_reads", "MT_reads"]))
    print("\t".join([
        args.sample,
        str(total_reads),
        str(mapped_reads),
        str(chrM_reads)
    ]))

if __name__ == "__main__":
    main()
