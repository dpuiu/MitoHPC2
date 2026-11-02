#!/usr/bin/env python3

import argparse
import re
import sys
import signal
from typing import TextIO, Optional, TypedDict


# Prevent BrokenPipeError when piping to head/tail
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


class Args(TypedDict):
    input: Optional[str]
    output: Optional[str]
    min_score: int
    add: int
    rmsuffix: bool
    offset: bool
    ed: bool


def parse_offset(chrom: str, start: int, end: int) -> tuple[str, int, int]:
    m = re.match(r"^(.+):(\d+)-(\d+)$", chrom)
    if m is not None:
        chrom = m.group(1)
        offset = int(m.group(2))
        start += offset
        end += offset
    return chrom, start, end


def process_bed(
    infile: TextIO,
    outfile: TextIO,
    add: int = 0,
    offset: bool = False,
    rmsuffix: bool = False,
    min_score: int = 0,
    ed: bool = False
) -> None:
    """
    Convert BED3–5 input to BED6 output with optional transformations:

    Features:
    - Supports BED lines with 3 to 5 columns
    - Optional coordinate offset parsing (chr:start-end)
    - Optional expansion of intervals (--add)
    - Optional removal of suffix from name (--rmsuffix)
    - Optional strand correction (--ed)
    - Filters by minimum BED score (--min)

    Input:
        BED from infile:
        chrom \t start \t end [\t name \t score]

    Output:
        BED6:
        chrom \t start \t end \t name \t score \t strand
    """

    for raw_line in infile:
        line: str = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        fields: list[str] = line.split("\t")

        if len(fields) < 3:
            raise ValueError(f"Invalid BED line: {line}")

        chrom: str = fields[0]
        start: int = int(fields[1])
        end: int = int(fields[2])

        if offset:
            chrom, start, end = parse_offset(chrom, start, end)

        if add != 0:
            start -= add
            end += add

        name: str = "."
        score: str = "."
        strand: str = "."

        if len(fields) >= 4:
            name = fields[3]
            if rmsuffix and "/" in name:
                name = name.split("/")[0]

        if len(fields) >= 5:
            score = fields[4]

        if len(fields) >= 6:
            strand = fields[5]

        if len(fields) == 3:
            name = f"{chrom}:{start}-{end}"

        if len(fields) == 4:
            score = str(end - start)

        if ed:
            try:
                score = str((end - start) - int(score))
            except Exception:
                pass

        if start > end:
            start, end = end, start
            strand = "-"

        try:
            if abs(int(score)) < min_score:
                continue
        except Exception:
            pass

        outfile.write(
            "\t".join([
                chrom,
                str(start),
                str(end),
                name,
                score,
                strand
            ]) + "\n"
        )


def parse_args() -> Args:
    parser = argparse.ArgumentParser(description="Convert BED3-5 to BED6")

    parser.add_argument("-i", "--input", type=str, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)

    parser.add_argument("--min", dest="min_score", type=int, default=0)
    parser.add_argument("--add", type=int, default=0)
    parser.add_argument("--rmsuffix", action="store_true")
    parser.add_argument("--offset", action="store_true")
    parser.add_argument("--ed", action="store_true")

    args = parser.parse_args()

    return {
        "input": args.input,
        "output": args.output,
        "min_score": args.min_score,
        "add": args.add,
        "rmsuffix": args.rmsuffix,
        "offset": args.offset,
        "ed": args.ed,
    }


def main() -> None:
    args: Args = parse_args()

    infile: TextIO
    outfile: TextIO

    infile = open(args["input"]) if args["input"] is not None else sys.stdin
    outfile = open(args["output"], "w") if args["output"] is not None else sys.stdout

    try:
        process_bed(
            infile=infile,
            outfile=outfile,
            add=args["add"],
            offset=args["offset"],
            rmsuffix=args["rmsuffix"],
            min_score=args["min_score"],
            ed=args["ed"]
        )
    finally:
        if args["input"] is not None:
            infile.close()
        if args["output"] is not None:
            outfile.close()


if __name__ == "__main__":
    main()
