#!/usr/bin/env python3
"""Simulate a small stranded paired-end RNA-seq test set with known UTRs.

Writes into the given directory:
  genome.fa, annotation.gff3 (CDS-only gene models, Liftoff-like IDs gene-X / rna-X)
  truth_utr_lengths.tsv (gene, utr5, utr3; IDs without the 'gene-' prefix)
  libA/libB_R{1,2}.fastq.gz (reverse-stranded, dUTP-like), libC_R{1,2}.fastq.gz (forward-stranded)
Features exercised: genes on both strands, a CDS intron, a 5' UTR intron, a close same-strand
tandem pair (neighbour cap), a tRNA, and a 'chrM' contig meant for --exclude_chroms.
Transcript ends are heterogeneous (3' ends spread over the last 30% of the UTR) so coverage
tapers, as in real data.
"""
import gzip, os, random, sys

rng = random.Random(11)
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
rand = lambda n: "".join(rng.choice("ACGT") for _ in range(n))


def orf(n_codons):
    stops = {"TAA", "TAG", "TGA"}
    cods = []
    while len(cods) < n_codons:
        c = rand(3)
        if c not in stops:
            cods.append(c)
    return "ATG" + "".join(cods) + "TAA"


chrA = list(rand(70000))
genes = []   # dict(name, strand, cds:[(s,e)], u5, u3, u5_intron:(s,e)|None)
x = 1000
specs = [  # (strand, u5, u3, cds_intron, u5_intron, gap_after)
    ("+", 60, 150, False, False, 600), ("+", 40, 220, True, False, 600), ("-", 80, 120, False, False, 600),
    ("+", 120, 180, False, True, 600), ("-", 30, 260, True, False, 600), ("+", 50, 140, False, False, 600),
    ("-", 100, 200, False, False, 600), ("+", 70, 160, False, False, 600), ("-", 45, 110, False, False, 600),
    ("+", 90, 240, False, False, 600), ("-", 60, 170, False, False, 600), ("+", 35, 130, False, False, 600),
    ("+", 55, 300, False, False, 40),   # tandem: next gene starts 40 bp after this 3' end -> capped
    ("+", 60, 150, False, False, 600), ("-", 75, 190, False, False, 600), ("+", 65, 150, False, False, 600),
]
for i, (sd, u5, u3, cint, uint, gap) in enumerate(specs):
    name = f"YT{i:03d}W" if sd == "+" else f"YT{i:03d}C"
    cds_seq = orf(rng.randint(150, 300))
    # transcript layout in transcript orientation: [u5 exon2]{intron}[u5 exon1? ] ... keep simple:
    u5_seq, u3_seq = rand(u5), rand(u3)
    intr5 = ("GT" + rand(146) + "AG") if uint else ""
    intrC = ("GT" + rand(96) + "AG") if cint else ""
    # 5' UTR intron splits the 5' UTR 20 bp upstream of the start codon
    if uint:
        tx_genomic = u5_seq[:-20] + intr5 + u5_seq[-20:]
    else:
        tx_genomic = u5_seq
    c1 = len(cds_seq) // 2
    cds_genomic = cds_seq[:c1] + intrC + cds_seq[c1:]
    tx_genomic += cds_genomic + u3_seq
    L = len(tx_genomic)
    s = x
    seq = tx_genomic if sd == "+" else rc(tx_genomic)
    chrA[s - 1:s - 1 + L] = list(seq)
    # coordinates (genomic, 1-based)
    off5 = len(u5_seq) + len(intr5)
    if sd == "+":
        cs = s + off5
        if cint:
            cds = [(cs, cs + c1 - 1), (cs + c1 + len(intrC), cs + len(cds_genomic) - 1)]
        else:
            cds = [(cs, cs + len(cds_seq) - 1)]
        uin = (s + len(u5_seq) - 20, s + len(u5_seq) - 20 + len(intr5) - 1) if uint else None
    else:
        e = s + L - 1
        ce = e - off5
        if cint:
            cds = [(ce - len(cds_genomic) + 1, ce - c1 - len(intrC)), (ce - c1 + 1, ce)]
        else:
            cds = [(ce - len(cds_seq) + 1, ce)]
        uin = None
    genes.append(dict(name=name, strand=sd, cds=cds, u5=u5, u3=u3, uintron=uin,
                      tx=u5_seq + cds_seq + u3_seq, start=s, end=s + L - 1))
    x = s + L + gap
# tRNA (no CDS) and a chrM gene
trna_s = x + 200
chrA[trna_s - 1:trna_s + 71] = list(rand(72))
chrM = list(rand(6000))
m_cds = orf(200)
chrM[1000:1000 + len(m_cds)] = list(m_cds)

with open(f"{out}/genome.fa", "w") as o:
    for nm, sq in (("chrA", chrA), ("chrM", chrM)):
        s = "".join(sq)
        o.write(f">{nm}\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")

with open(f"{out}/annotation.gff3", "w") as o, open(f"{out}/truth_utr_lengths.tsv", "w") as t:
    o.write("##gff-version 3\n")
    t.write("gene\tutr5\tutr3\n")
    for g in genes:
        n, sd = g["name"], g["strand"]
        lo, hi = g["cds"][0][0], g["cds"][-1][1]
        o.write(f"chrA\tLiftoff\tgene\t{lo}\t{hi}\t.\t{sd}\t.\tID=gene-{n};Name={n};gene_biotype=protein_coding\n")
        o.write(f"chrA\tLiftoff\tmRNA\t{lo}\t{hi}\t.\t{sd}\t.\tID=rna-{n};Parent=gene-{n}\n")
        for k, (a, b) in enumerate(g["cds"], 1):
            o.write(f"chrA\tLiftoff\texon\t{a}\t{b}\t.\t{sd}\t.\tID=exon-{n}-{k};Parent=rna-{n}\n")
        ph = 0
        segs = g["cds"] if sd == "+" else g["cds"][::-1]
        for a, b in segs:
            o.write(f"chrA\tLiftoff\tCDS\t{a}\t{b}\t.\t{sd}\t{ph}\tID=cds-{n};Parent=rna-{n}\n")
            ph = (3 - ((b - a + 1 - ph) % 3)) % 3
        t.write(f"{n}\t{g['u5']}\t{g['u3']}\n")
    o.write(f"chrA\tLiftoff\tgene\t{trna_s}\t{trna_s + 71}\t.\t+\t.\tID=gene-tT01;gene_biotype=tRNA\n")
    o.write(f"chrA\tLiftoff\ttRNA\t{trna_s}\t{trna_s + 71}\t.\t+\t.\tID=rna-tT01;Parent=gene-tT01\n")
    o.write(f"chrA\tLiftoff\texon\t{trna_s}\t{trna_s + 71}\t.\t+\t.\tID=exon-tT01-1;Parent=rna-tT01\n")
    me = 1000 + len(m_cds)
    o.write(f"chrM\tLiftoff\tgene\t1001\t{me}\t.\t+\t.\tID=gene-Q0001;gene_biotype=protein_coding\n")
    o.write(f"chrM\tLiftoff\tmRNA\t1001\t{me}\t.\t+\t.\tID=rna-Q0001;Parent=gene-Q0001\n")
    o.write(f"chrM\tLiftoff\texon\t1001\t{me}\t.\t+\t.\tID=exon-Q0001-1;Parent=rna-Q0001\n")
    o.write(f"chrM\tLiftoff\tCDS\t1001\t{me}\t.\t+\t0\tID=cds-Q0001;Parent=rna-Q0001\n")


def fragments(tx, u5, u3, n):
    L = len(tx)
    for _ in range(n):
        tss = rng.randint(0, max(0, int(u5 * 0.15)))            # 5' ends within the first 15% of the UTR
        pas = L - rng.randint(0, int(u3 * 0.3))                  # 3' ends spread over last 30% of the UTR
        mol = tx[tss:pas]
        fl = max(120, int(rng.gauss(200, 30)))
        if fl >= len(mol):
            continue
        a = rng.randint(0, len(mol) - fl)
        yield mol[a:a + fl]


def write_lib(name, stranded, nfrag):
    q = "I" * 100
    with gzip.open(f"{out}/{name}_R1.fastq.gz", "wt") as r1, gzip.open(f"{out}/{name}_R2.fastq.gz", "wt") as r2:
        k = 0
        for g in genes:
            for fr in fragments(g["tx"], g["u5"], g["u3"], nfrag):
                k += 1
                sense, anti = fr[:100], rc(fr[-100:])
                a, b = (anti, sense) if stranded == "reverse" else (sense, anti)
                r1.write(f"@{name}_{k}/1\n{a}\n+\n{q}\n"); r2.write(f"@{name}_{k}/2\n{b}\n+\n{q}\n")


for lib, st in (("libA", "reverse"), ("libB", "reverse"), ("libC", "forward")):
    write_lib(lib, st, 2500)
print(f"wrote {len(genes)} genes, 3 libraries to {out}")
