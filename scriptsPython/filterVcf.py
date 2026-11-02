#!/usr/bin/env python3

import argparse
import sys
import re
import signal
from typing import TextIO, Optional, TypedDict


# Prevent BrokenPipeError (e.g. piping to head)
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


class Args(TypedDict):
    input: Optional[str]
    output: Optional[str]
    percent: float
    depth: int
    sample: Optional[str]
    source: Optional[str]
    header: Optional[str]
    suspicious: Optional[str]


def load_suspicious(path: Optional[str]) -> set[str]:
    s: set[str] = set()
    if path:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    s.add(line.split("\t")[0])
    return s


def process_vcf(
    infile: TextIO,
    outfile: TextIO,
    percent: float = 0.0,
    depth: int = 0,
    sample: Optional[str] = None,
    header: Optional[str] = None,
    suspicious_path: Optional[str] = None
) -> None:

    suspicious = load_suspicious(suspicious_path)

    # --- custom header mode ---
    if header:
        with open(header) as f:
            for hdr_line in f:
                outfile.write(hdr_line)

        if sample:
            outfile.write(f"##sample={sample}\n")

        outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

    for raw_line in infile:
        line: str = raw_line.rstrip("\n")

        # header lines
        if line.startswith("#"):
            if not header:
                outfile.write(line + "\n")
            continue

        fields: list[str] = line.split("\t")

        if len(fields) < 10:
            raise ValueError(f"Invalid VCF line: {line}")

        chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample_field = fields[:10]

        # --- normalization check ---
        if "," in alt:
            raise ValueError(f"NORMALIZATION ERROR:{line}")

        # uppercase if percent > 0
        if percent > 0:
            ref = ref.upper()
            alt = alt.upper()

        # --- CASE 1: FORMAT == SM ---
        if fmt == "SM":

            # AF from INFO
            m = re.search(r"AF=([0-9]*\.?[0-9]+)", info)
            if m:
                af = float(m.group(1))
                if af < percent:
                    continue
                elif af > 1 - percent:
                    af = 1.0

                info = re.sub(r"AF=[0-9]*\.?[0-9]+", f"AF={af}", info)

            # DP filter
            m = re.search(r"DP=(\d+)", info)
            if m and int(m.group(1)) < depth:
                continue

            if sample_field in suspicious:
                continue

            outfile.write("\t".join([
                chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample_field
            ]) + "\n")

        # --- CASE 2: general FORMAT ---
        else:

            # add SM if missing
            if "SM=" not in info and sample:
                info = f"SM={sample}"
                if "INDEL" not in info and (alt == "*" or len(ref) != len(alt)):
                    info += ";INDEL"

            keys = fmt.split(":")
            vals = sample_field.split(":")
            h = dict(zip(keys, vals))

            # AF handling
            af_str = h.get("AF", "1")
            m = re.search(r"([0-9]*\.?[0-9]+)", af_str)
            af = float(m.group(1)) if m else 1.0

            if af > 1 - percent:
                af = 1.0
            elif af < percent:
                continue

            # GT logic
            if af == 1:
                gt = "1"
            elif af < 0.5:
                gt = "0/1"
            else:
                gt = "1/0"

            dp = int(h.get("DP", "0"))
            if dp < depth:
                continue

            ad = h.get("AD", ".")

            # suspicious via SM in INFO
            m = re.search(r"SM=([^;\s]+)", info)
            if m and m.group(1) in suspicious:
                continue

            fmt = "GT:DP:AD:AF"
            sample_field = f"{gt}:{dp}:{ad}:{af}"

            outfile.write("\t".join([
                chrom, pos, vid, ref, alt, qual, flt, info, fmt, sample_field
            ]) + "\n")


def parse_args() -> Args:
    parser = argparse.ArgumentParser(description="Filter VCF file")

    parser.add_argument("-i", "--input", type=str, default=None)
    parser.add_argument("-o", "--output", type=str, default=None)

    parser.add_argument("-percent", type=float, default=0.0)
    parser.add_argument("-depth", type=int, default=0)
    parser.add_argument("-sample", "--Run", dest="sample", type=str)
    parser.add_argument("-source", type=str)
    parser.add_argument("-header", type=str)
    parser.add_argument("-suspicious", type=str)

    args = parser.parse_args()

    if not (0 <= args.percent <= 1):
        raise ValueError("ERROR: percent must be >=0 and <=1")

    return {
        "input": args.input,
        "output": args.output,
        "percent": args.percent,
        "depth": args.depth,
        "sample": args.sample,
        "source": args.source,
        "header": args.header,
        "suspicious": args.suspicious,
    }


def main() -> None:
    args: Args = parse_args()

    infile: TextIO
    outfile: TextIO

    infile = open(args["input"]) if args["input"] else sys.stdin
    outfile = open(args["output"], "w") if args["output"] else sys.stdout

    try:
        process_vcf(
            infile=infile,
            outfile=outfile,
            percent=args["percent"],
            depth=args["depth"],
            sample=args["sample"],
            header=args["header"],
            suspicious_path=args["suspicious"]
        )
    finally:
        if args["input"]:
            infile.close()
        if args["output"]:
            outfile.close()


if __name__ == "__main__":
    main()
