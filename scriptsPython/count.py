#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def main():
    """
    Count values in a specified column of a tab-delimited file.

    Supports optional weighted counting and rounding.

    Input:
        tab-delimited text from stdin

    Example:
        cat file.tsv | python count.py -i 3 -j 4

    Options:
        -i  Column index used as key (default: 0)
        -j  Column index used as value (optional; default = count occurrences)
        --min   minimum count filter
        --max   maximum count filter
        --round round key values
    """

    parser = argparse.ArgumentParser(
        description="Count values in a specific column"
    )

    parser.add_argument("-i", type=int, default=0)
    parser.add_argument("-j", type=int)
    parser.add_argument("--min", type=int)
    parser.add_argument("--max", type=int)
    parser.add_argument("--round", type=int)

    args = parser.parse_args()

    count = defaultdict(int)

    for line in sys.stdin:
        if not line.strip() or line.startswith("#"):
            continue

        f = line.rstrip().split("\t")

        if len(f) <= args.i:
            continue
        if args.j is not None and len(f) <= args.j:
            continue

        key = f[args.i]

        # optional rounding of key
        if args.round:
            try:
                key = str(int(float(key) * args.round + 0.5) / args.round)
            except:
                pass

        if args.j is not None:
            try:
                count[key] += int(float(f[args.j]))
            except:
                continue
        else:
            count[key] += 1

    for k in sorted(count, key=lambda x: -count[x]):
        v = count[k]

        if args.min is not None and v < args.min:
            continue
        if args.max is not None and v > args.max:
            continue

        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
