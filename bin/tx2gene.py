#!/usr/bin/env python3
"""transcript -> gene table for tximport: transcript IDs from a FASTA (gffread -w output),
gene = the transcript's GFF3 Parent (parentless transcripts map to themselves).
Usage: tx2gene.py annotation.gff3 transcripts.fa > tx2gene.tsv"""
import re, sys
gff, fa = sys.argv[1], sys.argv[2]
tx = [l[1:].split()[0] for l in open(fa) if l.startswith(">")]
parent = {}
for l in open(gff):
    if l.startswith("##FASTA"):
        break
    if l.startswith("#") or l.count("\t") < 8:
        continue
    a = l.rstrip("\n").split("\t")[8]
    i = re.search(r"(?:^|;)ID=([^;]+)", a); p = re.search(r"(?:^|;)Parent=([^;,]+)", a)
    if i:
        parent[i.group(1)] = p.group(1) if p else None
miss = 0
for t in tx:
    if t not in parent:
        miss += 1
    print(f"{t}\t{parent.get(t) or t}")
print(f"tx2gene: {len(tx)} transcripts, {miss} not found as GFF IDs", file=sys.stderr)
