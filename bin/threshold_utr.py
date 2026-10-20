#!/usr/bin/env python3
"""Threshold-walk UTR caller, optionally calibrated on a table of known UTR lengths, with
optional utrrnaseq runs scored against it.

For each representative gene (from 02_make_utr_inputs.py):
  1. CDS coverage = median strand-specific pooled coverage over the CDS span.
  2. Build the outward coverage profile from the base next to the start codon (5') or stop
     codon (3'), in transcript coordinates: strong same-strand junctions met on the way are
     spliced over (UTR introns). The walk is capped at UTR_LIMIT transcript bp and at the
     next same-strand start/stop codon (same cap utrrnaseq uses).
  3. Smooth with a centred running mean (--smooth bp) and stop at the first position where the
     smoothed coverage < max(f * CDS coverage, --floor). UTR length = bp walked before that.
Calibration (only with --ref): genes with a reference length are split at random (fixed seed)
into calibration and validation halves; f is chosen separately for 5' and 3' by minimum median
|called - reference| on the calibration half (unless --f5/--f3 fix it). All reported comparisons
use the validation half, every method scored on the same genes (equal N); 95% CIs bootstrap
over genes. Without --ref, --f5/--f3 are used as given (defaults 0.05 / 0.03, calibrated on
S. cerevisiae Y55 against S288C TIF-seq, Pelechano et al. 2013).
--ref format: TSV with header; columns gene, utr5, utr3 (NA allowed). Gene IDs are matched
after removing --ref-id-strip (regex) from the annotation's gene IDs.

Outputs in --utrdir: UTRs.threshold.gff (utrrnaseq-style features), utr_calls.tsv (per gene,
all methods), calibration.tsv (grid scores). Stdlib only.
"""
import argparse, bisect, os, random, re, statistics as st
from array import array
from collections import defaultdict

GRID = [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.5]


def load_sizes(fa):
    sizes, name = {}, None
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]; sizes[name] = 0
            else:
                sizes[name] += len(line.strip())
    return sizes


def load_wig(path, sizes):
    cov, chrom = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith("track"):
                continue
            if line.startswith("variableStep"):
                chrom = re.search(r"chrom=(\S+)", line).group(1)
                cov[chrom] = array("i", bytes(4 * (sizes[chrom] + 2)))
                continue
            p, v = line.split()
            cov[chrom][int(p)] = int(float(v))
    return cov


def gid(attr, key="gene_id"):
    m = re.search(key + r' "([^"]+)"', attr)
    return m.group(1) if m else None


def load_utrrnaseq(prefix):
    """prefix e.g. <utrdir>/UTRs -> reads UTRs.plus.gff + UTRs.minus.gff; returns {gene: {'5': len, '3': len}}"""
    out = defaultdict(lambda: {"5": 0, "3": 0})
    for s in ("plus", "minus"):
        fn = f"{prefix}.{s}.gff"
        if not os.path.exists(fn):
            raise SystemExit(f"missing {fn}")
        with open(fn) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[2] in ("5'-UTR", "3'-UTR"):
                    out[gid(f[8])][f[2][0]] += int(f[4]) - int(f[3]) + 1
    return out


def outward_profile(cv, p0, d, cap, limit, introns, chrom, sd, min_jr, jfrac):
    """Positions/coverage walking outward from p0 in direction d (+1/-1), splicing over introns.
    Returns (positions, coverages, set of intron (s,e) used)."""
    pos, cov, used = [], [], []
    p = p0
    n = len(cv) - 2
    while len(pos) < limit and 1 <= p <= n and (p - cap) * d < 0:
        pos.append(p); cov.append(cv[p])
        nxt = p + d
        key = (chrom, sd, p, d)
        if key in introns:
            s, e, m = introns[key]
            inside = cv[s:s + 20] if d == 1 else cv[max(e - 19, 1):e + 1]
            if m >= min_jr and m >= jfrac * max(cv[p], 1) and (st.median(inside) if len(inside) else 0) < 0.5 * max(cv[p], 1):
                nxt = (e + 1) if d == 1 else (s - 1)
                used.append((s, e))
        p = nxt
    return pos, cov, used


def smooth(xs, w):
    if w <= 1 or not xs:
        return list(map(float, xs))
    h = w // 2
    pre = [0]
    for x in xs:
        pre.append(pre[-1] + x)
    out = []
    for i in range(len(xs)):
        a, b = max(0, i - h), min(len(xs), i + h + 1)
        out.append((pre[b] - pre[a]) / (b - a))
    return out


def call_len(sm, thr):
    for i, v in enumerate(sm):
        if v < thr:
            return i
    return len(sm)


def boot_ci(vals, fn, n=1000, seed=7):
    if len(vals) < 2:
        return ("NA", "NA")
    rng = random.Random(seed)
    stats = sorted(fn([vals[rng.randrange(len(vals))] for _ in vals]) for _ in range(n))
    return (round(stats[int(0.025 * n)], 1), round(stats[int(0.975 * n) - 1], 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--utrdir", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--ref", help="TSV of known UTR lengths: gene, utr5, utr3 (header row)")
    ap.add_argument("--ref-id-strip", default=r"^gene-", help="regex removed from GFF gene IDs before matching --ref")
    ap.add_argument("--limit", type=int, default=700)
    ap.add_argument("--smooth", type=int, default=11, help="running-mean window (bp)")
    ap.add_argument("--floor", type=float, default=5, help="absolute coverage floor")
    ap.add_argument("--min-cds-cov", type=float, default=10, help="genes below this get no call")
    ap.add_argument("--min-junction-reads", type=int, default=5)
    ap.add_argument("--junction-frac", type=float, default=0.2, help="junction reads / local depth to splice over")
    ap.add_argument("--utrrnaseq", action="append", default=[],
                    help="label=prefix of a utrrnaseq run, e.g. w150=<utrdir>/UTRs (repeatable)")
    ap.add_argument("--f5", type=float, help="fixed 5' fraction (no calibration); default 0.05 without --ref")
    ap.add_argument("--f3", type=float, help="fixed 3' fraction (no calibration); default 0.03 without --ref")
    ap.add_argument("--summary", default=None, help="also write the printed summary to this file")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    D = a.utrdir

    sizes = load_sizes(a.genome)
    cov = {s: load_wig(os.path.join(D, f"cov.{s}.wig"), sizes) for s in ("plus", "minus")}

    rep = {}
    with open(os.path.join(D, "representatives.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rep[f[0]] = dict(tid=f[1], chrom=f[2], strand=f[3], cs=int(f[4]), ce=int(f[5]))
    codons, allpos = defaultdict(dict), defaultdict(list)
    for s in ("plus", "minus"):
        with open(os.path.join(D, f"codons.{s}.gtf")) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                codons[gid(f[8])][f[2]] = (int(f[3]), int(f[4]))
                allpos[(f[0], f[6])] += [int(f[3]), int(f[4])]
    for k in allpos:
        allpos[k].sort()

    # junction lookup keyed by (chrom, strand, exon-side boundary base, walk direction)
    introns = {}
    with open(os.path.join(D, "introns.gff")) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            s, e, sd = int(f[3]), int(f[4]), f[6]
            m = int(re.search(r"mult=(\d+)", f[8]).group(1))
            for k, v in (((f[0], sd, s - 1, 1), (s, e, m)), ((f[0], sd, e + 1, -1), (s, e, m))):
                if k not in introns or introns[k][2] < m:
                    introns[k] = v

    tif = {}
    if a.ref:
        with open(a.ref) as fh:
            hdr = next(fh).rstrip("\n").split("\t")
            ix = {k: hdr.index(k) for k in ("gene", "utr5", "utr3")}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                tif[f[ix["gene"]]] = {e: None if f[ix["utr" + e]] in ("NA", "") else int(f[ix["utr" + e]]) for e in "53"}
    else:
        a.f5 = 0.05 if a.f5 is None else a.f5
        a.f3 = 0.03 if a.f3 is None else a.f3
        if a.utrrnaseq:
            raise SystemExit("--utrrnaseq comparison needs --ref")

    # build per-gene outward profiles once
    G = {}
    for g, r in rep.items():
        if "start_codon" not in codons[g] or "stop_codon" not in codons[g]:
            continue
        c, sd = r["chrom"], r["strand"]
        cv = cov["plus" if sd == "+" else "minus"].get(c)
        if cv is None:
            continue
        span = cv[r["cs"]:r["ce"] + 1]
        cds_cov = st.median(span) if len(span) else 0
        (ss, se), (ps, pe) = codons[g]["start_codon"], codons[g]["stop_codon"]
        pos = allpos[(c, sd)]
        if sd == "+":
            p5, d5 = ss - 1, -1; i = bisect.bisect_left(pos, ss); cap5 = pos[i - 1] if i > 0 else 0
            p3, d3 = pe + 1, 1; j = bisect.bisect_right(pos, pe); cap3 = pos[j] if j < len(pos) else sizes[c] + 1
        else:
            p5, d5 = se + 1, 1; j = bisect.bisect_right(pos, se); cap5 = pos[j] if j < len(pos) else sizes[c] + 1
            p3, d3 = ps - 1, -1; i = bisect.bisect_left(pos, ps); cap3 = pos[i - 1] if i > 0 else 0
        prof = {}
        for end, p0, d, cap in (("5", p5, d5, cap5), ("3", p3, d3, cap3)):
            P, C, U = outward_profile(cv, p0, d, cap, a.limit, introns, c, sd, a.min_junction_reads, a.junction_frac)
            prof[end] = dict(pos=P, sm=smooth(C, a.smooth), introns=U, capped=len(P) < a.limit)
        G[g] = dict(r, cds_cov=cds_cov, prof=prof, orf=re.sub(a.ref_id_strip, "", g))

    def called(g, end, f):
        x = G[g]
        if x["cds_cov"] < a.min_cds_cov:
            return None
        L = call_len(x["prof"][end]["sm"], max(f * x["cds_cov"], a.floor))
        return L if L > 0 else None

    # calibration / validation split (genes with S288C length for that end, and cds cov >= min)
    rng = random.Random(a.seed)
    split = {}
    for end in ("5", "3"):
        eligible = sorted(g for g, x in G.items()
                          if x["orf"] in tif and tif[x["orf"]][end] is not None and x["cds_cov"] >= a.min_cds_cov)
        rng.shuffle(eligible)
        h = len(eligible) // 2
        split[end] = (eligible[:h], eligible[h:])

    calib_rows, best = [], {}
    for end in ("5", "3"):
        cal = split[end][0]
        fixed = a.f5 if end == "5" else a.f3
        scores = []
        if not cal and fixed is None:
            raise SystemExit(f"--ref has no usable {end}' lengths matching annotation genes; check --ref-id-strip")
        for f in ([fixed] if fixed is not None else GRID):
            if not cal:
                break
            d = [abs((called(g, end, f) or 0) - tif[G[g]["orf"]][end]) for g in cal]
            scores.append((st.median(d), f))
            calib_rows.append((end, f, len(cal), st.median(d),
                               sum(called(g, end, f) is not None for g in cal)))
        best[end] = min(scores)[1] if scores else fixed

    with open(os.path.join(D, "calibration.tsv"), "w") as o:
        o.write("end\tf\tN_calibration\tmedian_abs_diff_bp\tN_called\n")
        for row in calib_rows:
            o.write("\t".join(map(str, row)) + "\n")

    utrr = {lab: load_utrrnaseq(pre) for lab, pre in (x.split("=", 1) for x in a.utrrnaseq)}

    # write GFF + per-gene table
    gff = open(os.path.join(D, "UTRs.threshold.gff"), "w")
    tab = open(os.path.join(D, "utr_calls.tsv"), "w")
    labs = sorted(utrr)
    tab.write("\t".join(["gene", "transcript", "chrom", "strand", "cds_cov",
                         "thr_utr5", "thr_utr3", "cap5", "cap3", "n_utr_introns5", "n_utr_introns3"]
                        + [f"{l}_utr{e}" for l in labs for e in "53"] + ["ref_utr5", "ref_utr3", "split5", "split3"]) + "\n")
    sp_lab = {e: {g: "cal" for g in split[e][0]} | {g: "val" for g in split[e][1]} for e in "53"}
    for g in sorted(G, key=lambda k: (G[k]["chrom"], G[k]["cs"])):
        x = G[g]
        row = [g, x["tid"], x["chrom"], x["strand"], round(x["cds_cov"], 1)]
        L, nint = {}, {"5": 0, "3": 0}
        for end in ("5", "3"):
            L[end] = called(g, end, best[end])
            if L[end]:
                P = x["prof"][end]["pos"][:L[end]]
                lo, hi = min(P), max(P)
                used = [(s, e) for s, e in x["prof"][end]["introns"] if lo < s and e < hi]
                # split the span into exon pieces around used introns
                pieces, start = [], lo
                for s, e in sorted(used):
                    pieces.append((start, s - 1)); start = e + 1
                pieces.append((start, hi))
                feat = "5'-UTR" if end == "5" else "3'-UTR"
                grp = f'gene_id "{g}"; transcript_id "{x["tid"]}";'
                for s, e in pieces:
                    gff.write(f"{x['chrom']}\tthresholdUTR\t{feat}\t{s}\t{e}\t.\t{x['strand']}\t.\t{grp}\n")
                nint[end] = len(used)
                for s, e in used:
                    gff.write(f"{x['chrom']}\tthresholdUTR\tintron\t{s}\t{e}\t.\t{x['strand']}\t.\t{grp}\n")
        cap = {e: (L[e] is not None and L[e] == len(x["prof"][e]["pos"]) and x["prof"][e]["capped"]) for e in "53"}
        row += [L["5"] or "NA", L["3"] or "NA", int(cap["5"]), int(cap["3"]), nint["5"], nint["3"]]
        for l in labs:
            row += [utrr[l][g]["5"] or "NA", utrr[l][g]["3"] or "NA"]
        t = tif.get(x["orf"], {"5": None, "3": None})
        row += [t["5"] if t["5"] is not None else "NA", t["3"] if t["3"] is not None else "NA",
                sp_lab["5"].get(g, "-"), sp_lab["3"].get(g, "-")]
        tab.write("\t".join(map(str, row)) + "\n")
    gff.close(); tab.close()

    # summary on validation half, same genes for every method
    lines = []
    def out(x=""):
        print(x); lines.append(x)
    out(f"genes evaluated: {len(G)}  (CDS cov >= {a.min_cds_cov}: {sum(x['cds_cov'] >= a.min_cds_cov for x in G.values())})")
    out(f"f used: 5'={best['5']}  3'={best['3']}  " + ("(calibrated; grid scores in calibration.tsv)" if a.ref and (a.f5 is None or a.f3 is None) else "(fixed)"))
    for e in ("5", "3"):
        if a.ref and (a.f5 if e == "5" else a.f3) is None and best[e] in (GRID[0], GRID[-1]):
            out(f"WARNING: calibrated {e}' f = {best[e]} is at the edge of the grid {GRID[0]}-{GRID[-1]}; the optimum may lie outside it")
    for end in ("5", "3"):
        if not split[end][1]:
            n_all = sum(1 for g in G if called(g, end, best[end]))
            out(f"\n== {end}' UTR: threshold calls {n_all}/{len(G)} (no reference lengths, no accuracy scoring)")
            continue
        val = split[end][1]
        out(f"\n== {end}' UTR, validation half: N = {len(val)} genes (same genes for every method)")
        out(f"{'method':<14}{'N':>6}{'called':>8}{'med_called':>11}{'med_ref':>10}{'med|diff|':>10}{'95% CI':>16}{'<=50bp':>8}")
        methods = [("threshold", lambda g: called(g, end, best[end]))]
        methods += [(l, (lambda l: lambda g: utrr[l][g][end] or None)(l)) for l in labs]
        for name, fn in methods:
            calls = {g: fn(g) for g in val}
            d = [abs((calls[g] or 0) - tif[G[g]["orf"]][end]) for g in val]   # uncalled scored as 0 bp
            got = [g for g in val if calls[g]]
            mc = st.median(calls[g] for g in got) if got else "NA"
            ms = st.median(tif[G[g]["orf"]][end] for g in val)
            ci = boot_ci(d, st.median)
            out(f"{name:<14}{len(val):>6}{len(got):>8}{mc:>11}{ms:>10}{st.median(d):>10}{str(ci):>16}{sum(x <= 50 for x in d):>8}")
        n_all = sum(1 for g in G if called(g, end, best[end]))
        n_cap = sum(1 for g in G if called(g, end, best[end]) and called(g, end, best[end]) == len(G[g]["prof"][end]["pos"]) and G[g]["prof"][end]["capped"])
        out(f"  all genes: threshold calls {n_all}/{len(G)}; ended at the neighbour-codon cap: {n_cap}")
    if a.ref:
        out("\nUncalled UTRs are scored as length 0 in |diff|. CIs: 1000 bootstrap resamples of genes.")
    if a.summary:
        open(a.summary, "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
