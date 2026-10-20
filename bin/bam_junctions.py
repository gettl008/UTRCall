#!/usr/bin/env python3
"""Count spliced junctions per transcript strand from SAM records on stdin.

Input: SAM text (no header) from `samtools view -q MAPQ -f 2 -F 0xF0C <bam>`, i.e. already
restricted to uniquely mapped, properly paired primary alignments.
Strand: reverse-stranded (dUTP) libraries put R2 on the transcript strand (R1 opposite);
forward-stranded libraries put R1 on the transcript strand.
Output (stdout, TSV, 1-based inclusive intron coordinates): chrom, start, end, strand, reads.
Each read crossing a junction counts once (mates crossing the same junction count twice).
"""
import argparse, re, sys
from collections import Counter

OP = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strandedness", choices=["reverse", "forward"], required=True)
    a = ap.parse_args()
    cnt = Counter()
    for line in sys.stdin:
        f = line.split("\t", 6)
        cig = f[5]
        if "N" not in cig:
            continue
        flag = int(f[1])
        rev, r1 = bool(flag & 16), bool(flag & 64)
        # transcript strand of the fragment
        same_as_read = (r1 and a.strandedness == "forward") or (not r1 and a.strandedness == "reverse")
        read_strand = "-" if rev else "+"
        strand = read_strand if same_as_read else ("+" if read_strand == "-" else "-")
        pos = int(f[3])
        for n, op in OP.findall(cig):
            n = int(n)
            if op == "N":
                cnt[(f[2], pos, pos + n - 1, strand)] += 1
            if op in "MDN=X":
                pos += n
    for (c, s, e, sd), n in sorted(cnt.items()):
        sys.stdout.write(f"{c}\t{s}\t{e}\t{sd}\t{n}\n")


if __name__ == "__main__":
    main()
