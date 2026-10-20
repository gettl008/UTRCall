#!/usr/bin/env python3
"""Write threshold-walk UTRs into the merged GFF3.

Rules (as agreed):
  * UTR ends that ran into the next same-strand codon (cap5/cap3 = 1 in utr_calls.tsv) are dropped
    unless --keep-capped.
  * UTRs are applied to every mRNA of the locus whose CDS span (min..max) equals the
    representative transcript's CDS span; other mRNAs are left unchanged.
  * For each modified end, exon sequence beyond the CDS on that end is replaced by the called
    UTR (adjacent piece merged into the terminal CDS exon, upstream/downstream pieces added as
    new exons when a UTR intron was spliced over). The stop codon may sit just outside the CDS
    (<= 3 bp gap); it is kept inside the exon.
  * five_prime_UTR / three_prime_UTR features of modified mRNAs are regenerated from exons - CDS;
    mRNA and gene coordinates are extended to their children.
  * mRNA gains utr5_len / utr3_len / utr_source attributes.
Self-checks after writing: every CDS lies inside an exon, exons of a transcript do not overlap.
Stdlib only.
"""
import argparse, os, re, sys
from collections import defaultdict, Counter
from urllib.parse import unquote

UTR_TYPES = {"five_prime_UTR", "three_prime_UTR"}
TX_TYPES = {"mRNA", "transcript"}


def attrs(col9):
    d = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--utrdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--source-tag", default="threshold_walk_RNAseq")
    ap.add_argument("--keep-capped", action="store_true", help="keep UTR ends stopped by the next same-strand codon (dropped by default)")
    a = ap.parse_args()
    rep = Counter()

    # ---- inputs from the UTR step
    calls = {}
    with open(os.path.join(a.utrdir, "utr_calls.tsv")) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            calls[r["gene"]] = r
    pieces = defaultdict(lambda: {"5": [], "3": []})
    with open(os.path.join(a.utrdir, "UTRs.threshold.gff")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] in ("5'-UTR", "3'-UTR"):
                g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                pieces[g][f[2][0]].append((int(f[3]), int(f[4])))

    # ---- parse GFF3 keeping original lines
    lines, feats = [], []
    with open(a.gff) as fh:
        for line in fh:
            lines.append(line)
    fasta_at = next((i for i, l in enumerate(lines) if l.startswith("##FASTA")), len(lines))
    by_id, children = {}, defaultdict(list)
    for i in range(fasta_at):
        l = lines[i]
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = attrs(f[8])
        ft = dict(i=i, f=f, at=at, id=at.get("ID"), parents=at.get("Parent", "").split(",") if at.get("Parent") else [])
        feats.append(ft)
        if ft["id"]:
            by_id[ft["id"]] = ft
        for p in ft["parents"]:
            children[p].append(ft)

    # ---- decide which transcripts get which ends
    plan = {}   # tid -> {"5": pieces or None, "3": pieces or None, gene}
    for g, r in calls.items():
        ends = {}
        for e in ("5", "3"):
            if r[f"thr_utr{e}"] == "NA":
                continue
            if r[f"cap{e}"] == "1" and not a.keep_capped:
                rep[f"dropped_neighbour_cap_{e}"] += 1
                continue
            if pieces[g][e]:
                ends[e] = sorted(pieces[g][e])
        if not ends:
            continue
        rtid = r["transcript"]
        rcds = [c for c in children.get(rtid, []) if c["f"][2] == "CDS"]
        if not rcds:
            rep["skip_rep_no_cds_in_gff"] += 1
            continue
        span = (min(int(c["f"][3]) for c in rcds), max(int(c["f"][4]) for c in rcds))
        gene_ft = by_id.get(g)
        txs = [c for c in children.get(g, []) if c["f"][2] in TX_TYPES] if gene_ft else []
        if not txs:
            txs = [by_id[rtid]] if rtid in by_id else []
        for t in txs:
            tc = [c for c in children.get(t["id"], []) if c["f"][2] == "CDS"]
            if not tc:
                continue
            tspan = (min(int(c["f"][3]) for c in tc), max(int(c["f"][4]) for c in tc))
            if tspan == span:
                plan[t["id"]] = dict(ends=ends, gene=g, strand=t["f"][6], cds=span)
                rep["transcripts_modified"] += 1
                if t["id"] != rtid:
                    rep["extra_isoforms_modified"] += 1
            else:
                rep["isoforms_left_unchanged_cds_differs"] += 1
        for e in ends:
            rep[f"genes_with_utr{e}"] += 1

    # ---- rebuild exons for planned transcripts
    new_children = {}   # tid -> list of (type, s, e, attrs-string) replacing exons + UTR features
    tx_bounds = {}
    for tid, p in plan.items():
        t = by_id[tid]
        sd = p["strand"]; cs, ce = p["cds"]
        ucs, uce = cs, ce   # UTR boundaries; widened past a stop codon that sits outside the CDS
        ex = [c for c in children[tid] if c["f"][2] == "exon"]
        if ex:
            exons = sorted([[int(c["f"][3]), int(c["f"][4]), c["f"][8]] for c in ex])
        else:
            rep["transcripts_without_exons_built_from_cds"] += 1
            exons = sorted([[int(c["f"][3]), int(c["f"][4]), None] for c in children[tid] if c["f"][2] == "CDS"])
        had_utr = any(s < cs or e > ce for s, e, _ in exons)
        side = {"5": "low" if sd == "+" else "high", "3": "high" if sd == "+" else "low"}
        for end, pcs in p["ends"].items():
            if side[end] == "low":
                adj = max(pcs, key=lambda x: x[1])
                if not (cs - 4 <= adj[1] <= cs - 1):
                    rep["skip_end_not_adjacent"] += 1; continue
                if any(s < cs for s, e, _ in exons):
                    rep["existing_utr_replaced"] += 1
                exons = [[max(s, cs), e, at] for s, e, at in exons if e >= cs]
                exons[0][0] = adj[0]
                ucs = adj[1] + 1
                exons += [[s, e, None] for s, e in pcs if (s, e) != adj]
            else:
                adj = min(pcs, key=lambda x: x[0])
                if not (ce + 1 <= adj[0] <= ce + 4):
                    rep["skip_end_not_adjacent"] += 1; continue
                if any(e > ce for s, e, _ in exons):
                    rep["existing_utr_replaced"] += 1
                exons = [[s, min(e, ce), at] for s, e, at in exons if s <= ce]
                exons[-1][1] = adj[1]
                uce = adj[0] - 1
                exons += [[s, e, None] for s, e in pcs if (s, e) != adj]
            exons.sort()
        rep["transcripts_had_utr_before"] += had_utr
        chrom, src = t["f"][0], t["f"][1]
        out, n = [], 0
        for s, e, at in exons:
            n += 1
            if at is None:
                at = f"ID={tid}-utrexon{n};Parent={tid};source={a.source_tag}"
            out.append((chrom, src, "exon", s, e, sd, at))
        # UTR features = exon parts outside CDS
        k = 0
        for s, e, _ in exons:
            for us, ue in ((s, min(e, ucs - 1)), (max(s, uce + 1), e)):
                if us > ue:
                    continue
                five = (ue < ucs) == (sd == "+")
                typ = "five_prime_UTR" if five else "three_prime_UTR"
                k += 1
                out.append((chrom, src, typ, us, ue, sd, f"ID={tid}-{typ}{k};Parent={tid};source={a.source_tag}"))
        new_children[tid] = out
        L5 = sum(e - s + 1 for c in out if c[2] == "five_prime_UTR" for s, e in [(c[3], c[4])])
        L3 = sum(e - s + 1 for c in out if c[2] == "three_prime_UTR" for s, e in [(c[3], c[4])])
        tx_bounds[tid] = (min(c[3] for c in out if c[2] == "exon"), max(c[4] for c in out if c[2] == "exon"), L5, L3)

    # gene bounds = union of their transcripts' (new or old) bounds
    gene_bounds = {}
    for tid, (s, e, _, _) in tx_bounds.items():
        g = by_id[tid]["parents"][0] if by_id[tid]["parents"] else None
        if g in by_id:
            gs, ge = gene_bounds.get(g, (int(by_id[g]["f"][3]), int(by_id[g]["f"][4])))
            gene_bounds[g] = (min(gs, s), max(ge, e))

    # ---- write
    skip_lines = set()
    for tid in new_children:
        for c in children[tid]:
            if c["f"][2] == "exon" or c["f"][2] in UTR_TYPES:
                skip_lines.add(c["i"])
    with open(a.out, "w") as o:
        for i, l in enumerate(lines):
            if i in skip_lines:
                continue
            if i >= fasta_at or l.startswith("#") or not l.strip():
                o.write(l); continue
            f = l.rstrip("\n").split("\t")
            at = attrs(f[8]) if len(f) >= 9 else {}
            fid = at.get("ID")
            if fid in tx_bounds:
                s, e, L5, L3 = tx_bounds[fid]
                f[3], f[4] = str(s), str(e)
                f[8] = f[8].rstrip(";") + f";utr5_len={L5};utr3_len={L3};utr_source={a.source_tag}"
                o.write("\t".join(f) + "\n")
                for c in new_children[fid]:
                    o.write("\t".join([c[0], c[1], c[2], str(c[3]), str(c[4]), ".", c[5], ".", c[6]]) + "\n")
                continue
            if fid in gene_bounds:
                f[3], f[4] = map(str, gene_bounds[fid])
                o.write("\t".join(f) + "\n"); continue
            o.write(l)

    # ---- self-checks on the written file
    ex, cds = defaultdict(list), defaultdict(list)
    with open(a.out) as fh:
        for l in fh:
            if l.startswith("##FASTA"):
                break
            if l.startswith("#") or not l.strip():
                continue
            f = l.rstrip("\n").split("\t")
            par = attrs(f[8]).get("Parent", "")
            if par in plan and f[2] == "exon":
                ex[par].append((int(f[3]), int(f[4])))
            elif par in plan and f[2] == "CDS":
                cds[par].append((int(f[3]), int(f[4])))
    bad_cds = sum(1 for t, cl in cds.items() for s, e in cl if not any(xs <= s and e <= xe for xs, xe in ex.get(t, [])))
    bad_ov = sum(1 for t, el in ex.items() for (a1, b1), (a2, b2) in zip(sorted(el), sorted(el)[1:]) if a2 <= b1)
    rep["check_CDS_outside_exon"] = bad_cds
    rep["check_overlapping_exons"] = bad_ov
    with open(a.out + ".report.txt", "w") as o:
        for k in sorted(rep):
            o.write(f"{k}\t{rep[k]}\n")
    sys.stdout.write(open(a.out + ".report.txt").read())
    if bad_cds or bad_ov:
        sys.exit("Self-check failed: see report")


if __name__ == "__main__":
    main()
