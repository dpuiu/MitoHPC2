#!/usr/bin/env python3

import sys
import argparse
import pysam


def main():
    """
    Print SAM lines whose QNAME appears in a provided list file.

    Input:
        sam_file   : SAM file (streamed)
        list_file  : file with query names (first column used)

    Output:
        filtered SAM (headers preserved)

    Note:
        This version streams SAM as text using pysam.AlignmentFile,
        but does NOT build alignment objects or modify BAM structure.
    """

    if len(sys.argv) < 3:
        sys.stderr.write("ERROR: insufficient number of arguments\n")
        sys.exit(1)

    sam_path = sys.argv[1]
    list_path = sys.argv[2]

    # -------------------------
    # load query names
    # -------------------------
    queries = set()

    with open(list_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            queries.add(line.split("\t")[0])

    if not queries:
        return

    # -------------------------
    # stream SAM as text (like Perl)
    # -------------------------
    with pysam.AlignmentFile(sam_path, "r") as sam:
        for line in sam.fetch(until_eof=True):
            # header lines are not returned by fetch(until_eof=True)
            # so we must handle SAM headers manually is NOT needed for BAM,
            # but for SAM we reopen raw text safely:

            break  # handled below

    # fallback: raw streaming (correct for SAM)
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                print(line, end="")
                continue

            fields = line.rstrip().split("\t")
            if not fields:
                continue

            if fields[0] in queries:
                print(line, end="")


if __name__ == "__main__":
    main()
