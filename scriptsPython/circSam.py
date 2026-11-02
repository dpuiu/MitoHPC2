#!/usr/bin/env python3

import sys
import argparse
import re

def load_lengths(fai_file):
    """
    Load reference lengths from .fai file
    """
    lengths = {}
    with open(fai_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = re.match(r"^(\S+)\s+(\d+)", line)
            if m:
                lengths[m.group(1)] = int(m.group(2))
    return lengths


def main():
    """
    Circularize SAM alignments by applying an offset and splitting
    alignments that cross reference boundaries.

    Example:
        cat I.sam | script.py --ref_len ref.fai --offset N > O.sam
    """

    parser = argparse.ArgumentParser(
        description="Circularize SAM alignments"
    )

    parser.add_argument("--ref_len", required=True, help="Reference .fai file")
    parser.add_argument("--offset", type=int, default=0)

    args = parser.parse_args()

    lengths = load_lengths(args.ref_len)

    for line in sys.stdin:
        line = line.rstrip("\n")

        # -----------------------------
        # Header handling
        # -----------------------------
        if line.startswith("@SQ"):
            fields = line.split("\t")
            m = re.search(r"SN:(\S+)", line)
            if m and m.group(1) in lengths:
                for i, f in enumerate(fields):
                    if f.startswith("LN:"):
                        fields[i] = f"LN:{lengths[m.group(1)]}"
            print("\t".join(fields))
            continue

        if line.startswith("@"):
            print(line)
            continue

        # -----------------------------
        # Alignment lines
        # -----------------------------
        F = line.split("\t")

        if len(F) < 6:
            continue

        rname = F[2]
        pos = int(F[3]) if F[3].isdigit() else 0
        cigar = F[5]
        rnext = F[6]
        pnext = int(F[7]) if F[7].isdigit() else 0

        # skip if no match in CIGAR
        if "M" not in cigar:
            print(line)
            continue

        # apply offset
        pos += args.offset
        pnext += args.offset

        # adjust PNEXT (mate position)
        if pnext:
            if rnext == "=" and rname in lengths and pnext > lengths[rname]:
                pnext %= lengths[rname]
            elif rnext in lengths and pnext > lengths[rnext]:
                pnext %= lengths[rnext]

        F[3] = str(pos)
        F[7] = str(pnext)

        # -----------------------------
        # Circular split logic
        # -----------------------------
        if rname in lengths and pos <= lengths[rname]:
            ref_len = lengths[rname]

            cigar1 = ""
            cigar2 = ""
            trim1 = 0
            trim2 = 0
            cur_pos = pos - 1
            c = cigar

            # soft/hard clip at start
            m = re.match(r"^(\d+)([SH])(.+)", c)
            if m:
                cigar1 += m.group(1) + m.group(2)
                trim2 += int(m.group(1))
                c = m.group(3)

            while True:
                m = re.match(r"^(\d+)(\w)(.*)", c)
                if not m:
                    break

                length = int(m.group(1))
                op = m.group(2)
                rest = m.group(3)

                if cur_pos + length <= ref_len:
                    cigar1 += f"{length}{op}"
                    if op != "D":
                        trim2 += length

                elif not trim1 and op == "D":
                    cigar1 += f"{length}{op}"

                elif not trim1:
                    trim1 = cur_pos + length - ref_len
                    keep1 = length - trim1

                    if keep1 > 0:
                        trim2 += keep1
                        cigar1 += f"{keep1}{op}"

                    cigar2 = f"{trim2}S{trim1}{op}"

                else:
                    cigar2 += f"{length}{op}"
                    if op != "D":
                        trim1 += length

                if op != "I":
                    cur_pos += length

                c = rest

            if trim1:
                cigar1 += f"{trim1}S"

            # merge trailing soft clips
            m = re.match(r"^(.+?)(\d+)S(\d+)S$", cigar1)
            if m:
                cigar1 = m.group(1) + str(int(m.group(2)) + int(m.group(3))) + "S"

            # print split read
            if "M" in cigar2:
                out = F.copy()
                out[1] = str(int(F[1]) + 2048)
                out[3] = "1"
                out[5] = cigar2
                print("\t".join(out))

            if "M" in cigar1:
                F[5] = cigar1
            else:
                F = []

        # wrap POS if needed
        if F and F[3].isdigit() and rname in lengths:
            if int(F[3]) > lengths[rname]:
                F[3] = str(int(F[3]) % lengths[rname])

        if F:
            print("\t".join(F))


if __name__ == "__main__":
    main()
