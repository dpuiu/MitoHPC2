#!/usr/bin/env python3

import sys

# Prevent BrokenPipeError when piping to head/tail
import signal
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

def main():
    try:
        for line in sys.stdin:

            # remove FILTER/FORMAT lines
            if line.startswith("##FILTER") or line.startswith("##FORMAT"):
                continue

            # keep header lines unchanged
            if line.startswith("##"):
                sys.stdout.write(line)
                continue

            # keep column header line
            if line.startswith("#CHROM"):
                sys.stdout.write("\t".join(line.rstrip().split("\t")[:8]) + "\n")
                continue

            # data lines → keep only columns 0–7
            cols = line.rstrip().split("\t")
            sys.stdout.write("\t".join(cols[:8]) + "\n")

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
