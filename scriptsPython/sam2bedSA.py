#!/usr/bin/env python3

import sys
import argparse
import pysam
import re
from typing import TypedDict, Optional, List, Tuple

import signal

# Prevent BrokenPipeError when piping to head/tail
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

# =============================
# types
# =============================
class SARecord(TypedDict):
    ref: str
    pos: int
    strand: str
    cigar: str


CIGAR_OPS: dict[str, int] = {
    "M": 0,
    "I": 1,
    "D": 2,
    "N": 3,
    "S": 4,
    "H": 5,
    "P": 6,
    "=": 7,
    "X": 8,
}


# =============================
# utils
# =============================
def cigar_to_tuples(cigar: str) -> List[Tuple[int, int]]:
    return [
        (CIGAR_OPS[op], int(n))
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


def add_mate_suffix(read: pysam.AlignedSegment) -> str:
    name: str = read.query_name or ""

    if read.is_read1:
        return f"{name}/1"
    if read.is_read2:
        return f"{name}/2"
    return name


def cigar_end(start: int, cigartuples: List[Tuple[int, int]]) -> int:
    end = start

    for op, length in cigartuples:
        if op in (0, 2, 3):  # M, D, N
            end += length

    return end


def parse_sa_tag(sa_string: str) -> Optional[SARecord]:
    parts = sa_string.split(";")[0].split(",")

    if len(parts) < 6:
        return None

    return SARecord(
        ref=str(parts[0]),
        pos=int(parts[1]) - 1,
        strand=str(parts[2]),
        cigar=str(parts[3]),
    )


# =============================
# core logic
# =============================
def process_read(read: pysam.AlignedSegment) -> None:
    """
    Process one read and emit primary + SA alignment.
    """

    if read.is_unmapped or read.is_secondary:
        return

    if not read.has_tag("SA"):
        return

    sa_raw: object = read.get_tag("SA")
    sa = parse_sa_tag(str(sa_raw))

    if sa is None:
        return

    qname: str = add_mate_suffix(read)
    strand: str = "-" if read.is_reverse else "+"

    start: int = read.reference_start
    cigartuples: List[Tuple[int, int]] = read.cigartuples or []
    end: int = cigar_end(start, cigartuples)
    score: int = end - start

    sa_start: int = sa["pos"]
    sa_cigartuples: List[Tuple[int, int]] = cigar_to_tuples(sa["cigar"])
    sa_end: int = cigar_end(sa_start, sa_cigartuples)
    sa_score: int = sa_end - sa_start

    primary = (
        read.reference_name or ".",
        start,
        end,
        qname,
        score,
        strand,
        read.cigarstring or "*",
    )

    secondary = (
        sa["ref"],
        sa_start,
        sa_end,
        qname,
        sa_score,
        sa["strand"],
        sa["cigar"],
    )

    if start < sa_start:
        print("\t".join(map(str, primary)) + "\t" + "\t".join(map(str, secondary)))
    else:
        print("\t".join(map(str, secondary)) + "\t" + "\t".join(map(str, primary)))


# =============================
# CLI
# =============================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Identify split alignments (SA tags)"
    )
    return parser.parse_args()


# =============================
# main
# =============================
def main() -> None:
    parse_args()

    with pysam.AlignmentFile(sys.stdin, "r") as sam:
        for read in sam:
            process_read(read)


if __name__ == "__main__":
    main()
