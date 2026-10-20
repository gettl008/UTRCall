#!/usr/bin/env python3
"""Build UTR-calling inputs from a GFF3 annotation, the genome FASTA and per-library junction
counts (bam_junctions.py output: chrom, intron_start, intron_end, strand, reads).

Outputs (in --outdir):
  codons.plus.gtf / codons.minus.gtf  start_codon + stop_codon per gene (one representative
                                      transcript per gene = longest CDS), split by strand
  introns.gff                         pooled STAR junctions, GT-AG only, utrrnaseq format
  representatives.tsv                 gene, transcript, chrom, strand, CDS span, codon flags
  prep_report.txt                     counts of everything kept / skipped

Stdlib only. Coordinates are 1-based inclusive throughout.
"""
import argparse, os, re, sys
from collections import Counter, defaultdict
from urllib.parse import unquote

STOPS_NUC = {"TAA", "TAG", "TGA"}
STOPS_ANY = STOPS_NUC  # mito TGA=Trp, but a TGA "stop" there only shifts the anchor by 3 bp
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def attrs(col9):
    d = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = unquote(v.strip())
    return d


def parse_gff(path):
    """Return transcripts {tid: dict(gene, chrom, strand, cds=[(s,e)])} and gene types."""
    parent_of, ftype, cds = {}, {}, defaultdict(list)
    chrom_strand = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            fid = a.get("ID")
            par = a.get("Parent", "").split(",")[0] or None
            if f[2] == "CDS":
                if par:
                    cds[par].append((int(f[3]), int(f[4])))
                    chrom_strand[par] = (f[0], f[6])
                continue
            if fid:
                ftype[fid] = f[2]
                if par:
                    parent_of[fid] = par
    tx = {}
    for tid, segs in cds.items():
        gene = parent_of.get(tid, tid)
        tx[tid] = dict(gene=gene, chrom=chrom_strand[tid][0], strand=chrom_strand[tid][1],
                       cds=sorted(segs), gtype=ftype.get(gene, "?"), ttype=ftype.get(tid, "?"))
    return tx


def codon_at(seq, s, e, strand):
    sub = seq[s - 1:e]
    return rc(sub) if strand == "-" else sub


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--junctions", nargs="+", required=True, help="bam_junctions.py TSVs (one per library)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-intron-reads", type=int, default=5, help="min pooled unique reads per junction")
    ap.add_argument("--exclude-chroms", default="", help="comma-separated seqids to skip (e.g. mito)")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    rep = Counter()
    excl = {c for c in a.exclude_chroms.split(",") if c}
    genome = read_fasta(a.genome)
    tx = parse_gff(a.gff)
    rep["transcripts_with_CDS"] = len(tx)

    # one representative per gene: longest total CDS
    best = {}
    for tid, t in tx.items():
        if t["chrom"] in excl:
            rep["skip_excluded_chrom"] += 1
            continue
        if t["chrom"] not in genome:
            rep["skip_chrom_not_in_fasta"] += 1
            continue
        L = sum(e - s + 1 for s, e in t["cds"])
        if t["gene"] not in best or L > best[t["gene"]][0]:
            best[t["gene"]] = (L, tid)
    rep["genes_represented"] = len(best)

    out = {"+": open(os.path.join(a.outdir, "codons.plus.gtf"), "w"),
           "-": open(os.path.join(a.outdir, "codons.minus.gtf"), "w")}
    reptab = open(os.path.join(a.outdir, "representatives.tsv"), "w")
    reptab.write("gene\ttranscript\tchrom\tstrand\tcds_start\tcds_end\tstart_codon\tstop_codon\tstop_source\tgene_type\n")
    per_chrom = Counter()

    for gene, (L, tid) in sorted(best.items(), key=lambda x: (tx[x[1][1]]["chrom"], tx[x[1][1]]["cds"][0][0])):
        t = tx[tid]
        chrom, strand, segs = t["chrom"], t["strand"], t["cds"]
        if strand not in "+-":
            rep["skip_no_strand"] += 1
            continue
        seq = genome[chrom]
        cs, ce = segs[0][0], segs[-1][1]
        # start codon = first 3 nt of CDS in transcript orientation; stop = last 3 (or next 3)
        if strand == "+":
            first, last = segs[0], segs[-1]
            st = (cs, cs + 2)
            sp_in = (ce - 2, ce)
            sp_out = (ce + 1, ce + 3)
        else:
            first, last = segs[-1], segs[0]
            st = (ce - 2, ce)
            sp_in = (cs, cs + 2)
            sp_out = (cs - 3, cs - 1)
        if first[1] - first[0] < 2 or last[1] - last[0] < 2:
            rep["skip_codon_split_by_intron"] += 1
            continue
        sc = codon_at(seq, *st, strand)
        stop_in = codon_at(seq, *sp_in, strand)
        stop_out = codon_at(seq, *sp_out, strand) if sp_out[0] >= 1 and sp_out[1] <= len(seq) else ""
        if stop_in in STOPS_ANY:
            sp, src = sp_in, "in_CDS"
        elif stop_out in STOPS_ANY:
            sp, src = sp_out, "after_CDS"
        else:
            sp, src = sp_in, "none_used_CDS_end"
        rep["start_ATG" if sc == "ATG" else "start_nonATG"] += 1
        rep["stop_" + src] += 1
        stopc = codon_at(seq, *sp, strand)
        grp = f'gene_id "{gene}"; transcript_id "{tid}";'
        o = out[strand]
        o.write(f"{chrom}\tY55\tstart_codon\t{st[0]}\t{st[1]}\t0\t{strand}\t.\t{grp}\n")
        o.write(f"{chrom}\tY55\tstop_codon\t{sp[0]}\t{sp[1]}\t0\t{strand}\t.\t{grp}\n")
        reptab.write(f"{gene}\t{tid}\t{chrom}\t{strand}\t{cs}\t{ce}\t{sc}\t{stopc}\t{src}\t{t['gtype']}\n")
        per_chrom[chrom] += 1
        rep["genes_written_" + strand] += 1
    for o in out.values():
        o.close()
    reptab.close()

    # introns: pool per-library junction counts; keep GT-AG on the transcript strand
    sj = Counter()
    rep["junction_files"] = len(a.junctions)
    for fn in a.junctions:
        with open(fn) as fh:
            for line in fh:
                c, s, e, sd, n = line.rstrip("\n").split("\t")
                sj[(c, int(s), int(e), sd)] += int(n)
    motif = Counter()
    with open(os.path.join(a.outdir, "introns.gff"), "w") as fh:
        for (c, s, e, sd), n in sorted(sj.items()):
            if c not in genome or c in excl:
                continue
            seq = genome[c]
            m = seq[s - 1:s + 1] + seq[e - 2:e]
            if sd == "-":
                m = rc(seq[e - 2:e]) + rc(seq[s - 1:s + 1])
            motif[m] += 1
            if m != "GTAG" or n < a.min_intron_reads:
                continue
            fh.write(f"{c}\tRNAseq\tintron\t{s}\t{e}\t0\t{sd}\t.\tmult={n}_GTAG\n")
            rep["introns_written"] += 1
    rep["junctions_seen"] = len(sj)
    rep["junctions_GTAG"] = motif["GTAG"]

    with open(os.path.join(a.outdir, "prep_report.txt"), "w") as fh:
        for k, v in rep.items():
            fh.write(f"{k}\t{v}\n")
        fh.write("\n# genes per seqid\n")
        for c in sorted(per_chrom):
            fh.write(f"{c}\t{per_chrom[c]}\t{len(genome[c])}\n")
    sys.stdout.write(open(os.path.join(a.outdir, "prep_report.txt")).read())


if __name__ == "__main__":
    main()
