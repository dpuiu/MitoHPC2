#!/usr/bin/env python3

import sys
import pysam
import signal

# Prevent BrokenPipeError when piping to head/tail
signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def main():
    """
    Add AF, AD, GT FORMAT fields to bcftools VCF output using pysam.

    Functionality:
    - Skips REF alleles containing 'N'
    - Extracts DP and DP4-like values from INFO (if present)
    - Computes:
        AD = sum of ALT depths from DP4
        AF = allele frequency
    - Rewrites FORMAT as:
        GT:DP:AD:AF
    - Adds missing FORMAT definitions to header

    Input:
        VCF from stdin

    Output:
        Modified VCF to stdout
    """

    try:
        vcf_in = pysam.VariantFile(sys.stdin)

        header = vcf_in.header.copy()

        # Add FORMAT definitions (like Perl injection)
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.add_line('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">')
        header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">')

        vcf_out = pysam.VariantFile(sys.stdout, "w", header=header)

        for rec in vcf_in:

            # skip REF containing N
            if rec.ref and "N" in rec.ref:
                continue

            # -------------------------
            # INFO extraction
            # -------------------------
            dp = rec.info.get("DP")

            dp4 = rec.info.get("DP4")
            if not dp or not dp4:
                continue

            # DP4 = (ref_fw, ref_rev, alt_fw, alt_rev)
            try:
                r1, r2, a1, a2 = dp4
            except Exception:
                continue

            ad = a1 + a2
            total = r1 + r2 + a1 + a2

            af = round(ad / total, 3) if total > 0 else 0.0

            # -------------------------
            # GT extraction
            # -------------------------
            sample = list(rec.samples.values())[0]

            gt = sample.get("GT", (None,))
            if isinstance(gt, tuple):
                gt = "/".join("." if g is None else str(g) for g in gt)

            # -------------------------
            # rewrite FORMAT
            # -------------------------
            sample.clear()

            sample["GT"] = gt
            sample["DP"] = dp
            sample["AD"] = ad
            sample["AF"] = af

            rec.format = ["GT", "DP", "AD", "AF"]

            vcf_out.write(rec)

        vcf_out.close()

    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
