# utrcall

Adds RNA-seq-supported 5′ and 3′ UTRs to an existing gene annotation (for example a Liftoff
transfer that has CDS-only models). Built for *S. cerevisiae* Y55 and, theoretically, usable on any genome with
stranded paired-end RNA-seq data.

## Method

1. **Alignment** (FASTQ rows only): STAR genome index, then 2-pass alignment. Pre-aligned,
   coordinate-sorted BAMs are used as given.
2. **Junctions**: spliced reads from uniquely mapped, properly paired primary alignments
   (`--min_mapq`). The transcript strand comes from the library strandedness. Junctions are
   pooled across libraries, and GT–AG junctions with ≥ `--min_intron_reads` reads are kept.
3. **Coverage**: strand-specific depth pooled across all libraries. Libraries of each
   strandedness are pooled separately and then summed. Overlapping mates are counted once
   (`samtools depth -s`), and each seqid is processed separately to bound memory.
4. **Representative transcript per gene**: the longest CDS. The start and stop codons are taken
   from the genome sequence, and the stop codon can sit inside or just after the CDS.
5. **Threshold walk**: from each codon, the caller walks outward over smoothed coverage
   (`--smooth`) until it drops below `max(f × median CDS coverage, --floor)`.
   - It splices over strong same-strand UTR introns.
   - It stops at `--utr_limit` or at the next same-strand start/stop codon.
   - Genes below `--min_cds_cov` get no call.
6. **f**: either calibrated or fixed.
   - Calibrated: with `--ref_utr` (known UTR lengths), f is chosen separately for 5′ and 3′
     on a random half of the reference genes. All methods are then scored on the other half
     (equal N, bootstrap 95% CI over genes).
   - Fixed: without a reference, `--f5`/`--f3` are used as given. The defaults (0.05 / 0.03)
     were calibrated on Y55 against S288C TIF-seq (Pelechano et al. 2013).
7. **Merge into the GFF3**:
   - UTRs go on every mRNA of the locus whose CDS span equals the representative's.
   - UTR ends that ran into the next same-strand gene are dropped unless `--keep_capped`.
   - Exons, mRNA and gene coordinates are extended, and `five_prime_UTR`/`three_prime_UTR`
     features are regenerated.
   - Self-checks (CDS inside exons, no overlapping exons) and `gffread -E` must pass.

Optional modules:

- `--run_utrrnaseq`: AUGUSTUS `utrrnaseq` with smoothing windows of 150 and 50 bp, scored
  against the threshold caller. Needs `--ref_utr`.
- `--run_salmon`: Salmon 2.7.0 decoy-aware index on the UTR-updated annotation, a quant run per
  FASTQ library, and `tx2gene.tsv` for tximport.

On Y55 (42 libraries, held-out half of 2296 genes vs S288C TIF-seq), median |error| was
13 bp (5′) and 18 bp (3′) for the threshold walk, against 29 / 98 bp for utrrnaseq at
its defaults.

## Inputs

`--samplesheet` is a CSV with columns `sample,fastq_1,fastq_2,bam,strandedness`.

- Each row gives either a FASTQ pair or a BAM.
- `strandedness` is `reverse` (dUTP / TruSeq stranded, Salmon `ISR`) or `forward`.
  Unstranded libraries are rejected, because the method needs strand-specific coverage.
- Data must be paired-end.
- Relative paths are resolved against the samplesheet's directory.

See `assets/samplesheet_example.csv`.

`--genome` is the FASTA. `--gff` is the GFF3 with gene → mRNA → exon/CDS.

`--ref_utr` (optional) is a TSV with header `gene utr5 utr3` (NA allowed). Gene IDs are
matched after stripping `--ref_id_strip` (regex, default `^gene-`) from the GFF gene IDs.
`assets/s288c_tifseq_ypd_utr_lengths.tsv` holds S288C TIF-seq YPD lengths for 4818 ORFs,
derived from the SGD-hosted Pelechano 2013 most-abundant full-ORF transcripts.

## Running

Tested with Nextflow 25.04.6 and 26.04.6 (strict syntax; `nextflow lint` clean).

```bash
# local test (simulated 70 kb genome, 16 genes, 2 reverse FASTQ libraries + 1 forward BAM)
nextflow run main.nf -profile conda,test,local

# RDAC (SLURM): head job submits each task as its own job
sbatch --export=ALL,UTRCALL_DIR=$HOME/.local/soft/utrcall run_rdac.sbatch params/y55.yaml
sbatch --export=ALL,UTRCALL_DIR=$HOME/.local/soft/utrcall run_rdac.sbatch params/y55.yaml -resume
```

Profiles:

| Profile | What it sets |
|---|---|
| `conda` | Per-process envs from `envs/*.yml`, or existing envs via `--env_core` / `--env_augustus` / `--env_salmon` |
| `rdac` | SLURM, `-A rdac_acct`, `--exclude=node002`, queueSize 50 |
| `local` | Local executor |
| `test` | The bundled test data |

Key parameters are listed in `nextflow.config` and in `nextflow run main.nf --help`.

| Param | Default | Note |
|---|---|---|
| `min_mapq` | 255 | STAR unique; use 60 for HISAT2 BAMs |
| `exclude_chroms` | '' | e.g. the mitochondrion (polycistronic transcripts) |
| `utr_limit` | 700 | max UTR length walked |
| `intron_max` | 2000 | STAR intron/mate-gap limit; raise for larger genomes |
| `f5` / `f3` | 0.05 / 0.03 | used without `ref_utr`, or with `--fix_f` |
| `keep_capped` | false | keep UTR ends stopped by a neighbour |

## Outputs (`--outdir`)

```
<out_prefix>.utr.gff3                  annotation with UTRs
<out_prefix>.utr.gff3.report.txt       merge counts + self-checks
<out_prefix>.utr.gffread_check.txt     gffread -E log
utr/utr_summary.txt                    f used, call counts, accuracy vs reference (if given)
utr/utr_calls.tsv                      per gene: CDS coverage, UTR lengths, neighbour-cap flags,
                                       utrrnaseq lengths, reference lengths, calibration split
utr/calibration.tsv                    grid of f vs median |error| (with ref_utr)
utr/UTRs.threshold.gff                 raw UTR calls
utr/prep/                              representative codons, introns, prep_report.txt
utr/coverage/cov.{plus,minus}.wig      pooled strand-specific coverage
utr/utrrnaseq/                         (optional) utrrnaseq calls
star/                                  (FASTQ rows) BAMs, Log.final.out
salmon/                                (optional) quant/<sample>, tx2gene.tsv, transcripts.fa
pipeline_info/                         Nextflow report + trace
```

## Caveats

- UTR ends from short-read coverage are approximate. Yeast TSS and polyadenylation sites are
  heterogeneous, and the threshold picks the extent of the pooled population of transcripts,
  not a single isoform.
- Calibration against another strain or condition (S288C YPD for Y55) absorbs real
  differences, which is why scoring uses a held-out half.
- If calibrated f lands at the edge of the grid (0.005–0.5), `utr_summary.txt` warns. That
  happens with sharp-edged simulated data and is not expected on real libraries.
- Genes are anchored on the CDS of their representative transcript. Other isoforms of a locus
  with a different CDS span are left unchanged.
  
*Built with Claude Opus 5.5*
