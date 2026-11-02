#!/usr/bin/env python3

import sys
import argparse
import statistics


def main():
    """
    Generate summary statistics from numeric input.

    Input:
        One numeric value per line (int or float)

    Output:
        count, min, q1, median, q3, max, mean, stdev

    Example:
        cut -f3 I.cvg | python script.py > O.cvg.stat
    """

    parser = argparse.ArgumentParser(
        description="Generate summary statistics"
    )

    parser.add_argument("-sample", "--sample", help="Sample name (Run column)")
    args = parser.parse_args()

    values = []

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            values.append(float(line))
        except ValueError:
            sys.stderr.write(f"ERROR: {line}\n")
            sys.exit(1)

    if len(values) < 2:
        sys.exit(0)

    values.sort()

    # core stats
    n = len(values)
    min_v = values[0]
    max_v = values[-1]
    mean_v = statistics.mean(values)

    # stdev (sample standard deviation, like R's sd())
    stdev_v = statistics.stdev(values)

    # quartiles
    q1, median, q3 = statistics.quantiles(values, n=4, method="inclusive")

    # header
    if args.sample:
        print("Run\t", end="")
    print("count\tmin\tq1\tmedian\tq3\tmax\tmean\tstdev")

    # output
    if args.sample:
        print(f"{args.sample}\t", end="")

    print("\t".join([
        str(n),
        str(min_v),
        str(q1),
        str(median),
        str(q3),
        str(max_v),
        str(round(mean_v, 3)),
        str(round(stdev_v, 3))
    ]))


if __name__ == "__main__":
    main()
