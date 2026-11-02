#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict


def main():
    """
    Filter lines allowing up to MAX occurrences of unique values
    in a specified column.

    Behavior:
    - Keeps at most `max` occurrences per key
    - Preserves order (streaming)
    - Handles full-line mode if -i is not set

    Example:
        cat I.tab | python script.py -i 3 -max 2
    """

    parser = argparse.ArgumentParser(
        description="Limit occurrences of unique values per column"
    )

    parser.add_argument("-i", type=int, help="Column index (0-based)")
    parser.add_argument("-max", "--max", type=int, default=1)

    args = parser.parse_args()

    counts = defaultdict(int)

    for line in sys.stdin:
        line = line.rstrip("\n")

        if line == "":
            print("")
            continue

        # -------------------------
        # column mode
        # -------------------------
        if args.i is not None:
            if line.startswith("#") or line.startswith("@"):
                print(line)
                continue

            fields = line.split("\t")

            if args.i >= len(fields):
                continue

            key = fields[args.i]

        # -------------------------
        # full-line mode
        # -------------------------
        else:
            key = line

        counts[key] += 1

        if counts[key] <= args.max:
            print(line)


if __name__ == "__main__":
    main()
