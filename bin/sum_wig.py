#!/usr/bin/env python3
"""Sum variableStep wiggle files (same strand, different library groups) position-wise.
Usage: sum_wig.py genome.fa.fai out.wig in1.wig [in2.wig ...]"""
import re, sys
from array import array

fai, out, ins = sys.argv[1], sys.argv[2], sys.argv[3:]
sizes = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
order = list(sizes)
cov = {}
for fn in ins:
    chrom = None
    with open(fn) as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                chrom = re.search(r"chrom=(\S+)", line).group(1)
                if chrom not in cov:
                    cov[chrom] = array("l", bytes(8 * (sizes[chrom] + 2)))
                continue
            p, v = line.split()
            cov[chrom][int(p)] += int(float(v))
with open(out, "w") as o:
    name = re.sub(r"\.wig$", "", out.split("/")[-1])
    o.write(f"track type=wiggle_0 name={name}\n")
    for c in order:
        if c not in cov:
            continue
        o.write(f"variableStep chrom={c}\n")
        a = cov[c]
        o.write("".join(f"{i} {a[i]}\n" for i in range(1, sizes[c] + 1) if a[i]))
