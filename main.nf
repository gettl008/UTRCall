#!/usr/bin/env nextflow
/*
 * utrcall: add RNA-seq-supported 5'/3' UTRs to an existing gene annotation.
 *
 *   samplesheet (FASTQ pairs and/or coordinate-sorted BAMs, stranded PE)
 *     -> [STAR index + 2-pass align, FASTQ rows only]
 *     -> per-library junction counts (transcript strand from library strandedness)
 *     -> pooled strand-specific coverage (per strandedness group, summed)
 *     -> representative transcript per gene, start/stop codons, GT-AG introns
 *     -> threshold-walk UTR calls (optionally calibrated on known UTR lengths)
 *     -> UTRs merged into the GFF3 (+ gffread check)
 *   optional: utrrnaseq comparison runs; Salmon index + quant (FASTQ rows)
 */
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    utrcall v${workflow.manifest.version}

    Required:
      --samplesheet   CSV: sample,fastq_1,fastq_2,bam,strandedness  (fastq pair OR bam per row;
                      strandedness = reverse (dUTP/TruSeq stranded) | forward; relative paths are
                      resolved against the samplesheet's directory)
      --genome        genome FASTA
      --gff           GFF3 annotation (gene > mRNA > exon/CDS)
      --outdir        output directory

    Common options (defaults in nextflow.config):
      --exclude_chroms     comma-separated seqids skipped for UTR calling (e.g. mitochondrion)
      --min_mapq           unique-alignment MAPQ (STAR 255; HISAT2 60)          [${params.min_mapq}]
      --utr_limit          max UTR length walked (transcript bp)                [${params.utr_limit}]
      --ref_utr            TSV gene,utr5,utr3 of known UTR lengths -> calibrate f and score
      --ref_id_strip       regex removed from GFF gene IDs to match --ref_utr   [${params.ref_id_strip}]
      --f5 / --f3          fixed coverage-drop fractions (used when no --ref_utr) [${params.f5} / ${params.f3}]
      --keep_capped        keep UTR ends that ran into the next same-strand gene [${params.keep_capped}]
      --run_utrrnaseq      also run AUGUSTUS utrrnaseq for comparison (needs --ref_utr) [${params.run_utrrnaseq}]
      --run_salmon         Salmon 2.7 index + quant on FASTQ rows                [${params.run_salmon}]
      --star_index         prebuilt STAR index directory (skips index build)

    Profiles: conda, rdac (SLURM on DU RDAC), test
    """.stripIndent()
}
// ---------------------------------------------------------------- input parsing
// Strict-syntax friendly (Nextflow >= 25.10 parses scripts with the strict parser): all logic
// lives in functions or the workflow body, no top-level statements.

def refUtr() {
    // '--ref_utr null|false|""' on the command line means "no reference"
    def r = params.ref_utr
    return (r && !(r.toString() in ['null', 'false'])) ? r.toString() : null
}

def resolvePath(String p, sheetDir) {
    if (!p) return p
    return (p.startsWith('/') || p ==~ /^[a-z]+:\/\/.*/) ? p : sheetDir.resolve(p).toString()
}

def readSamplesheet(String path) {
    def sheet = file(path, checkIfExists: true)
    def sheetDir = sheet.toAbsolutePath().parent
    def rows = []
    def hdr = null
    def n = 0
    sheet.readLines().each { String line ->
        n += 1
        if (!line.trim() || line.startsWith('#')) return
        def f = line.split(',', -1).collect { String x -> x.trim() }
        if (hdr == null) { hdr = f; return }
        def r = [:]
        hdr.eachWithIndex { String h, int k -> r[h] = k < f.size() ? f[k] : '' }
        ['fastq_1', 'fastq_2', 'bam'].each { String k -> r[k] = resolvePath((r[k] ?: '') as String, sheetDir) }
        def fq = r.fastq_1 && r.fastq_2
        def bm = r.bam as boolean
        if (!r.sample) error "samplesheet line ${n}: empty sample"
        if (fq == bm) error "samplesheet line ${n} (${r.sample}): give either fastq_1+fastq_2 or bam, not both/neither"
        if (!(r.strandedness in ['reverse', 'forward'])) error "samplesheet line ${n} (${r.sample}): strandedness must be reverse or forward (unstranded libraries cannot give strand-specific coverage)"
        rows << r
    }
    if (!rows) error "samplesheet has no rows"
    if (rows.collect { r -> r.sample }.toSet().size() != rows.size()) error "duplicate sample names in samplesheet"
    return rows
}

def checkParams() {
    ['samplesheet', 'genome', 'gff', 'outdir'].each { String p ->
        if (!params[p]) { helpMessage(); error "Missing required --${p}" }
    }
    if (params.run_utrrnaseq && !refUtr()) error "--run_utrrnaseq needs --ref_utr (its comparison is scored against known UTR lengths)"
}

// ---------------------------------------------------------------- processes
process FAIDX {
    label 'small'
    input:  path genome
    output: tuple path(genome), path("${genome}.fai")
    script: "samtools faidx ${genome}"
}

process STAR_INDEX {
    label 'medium'
    publishDir "${params.outdir}/star", mode: 'copy', enabled: params.publish_star_index
    input:  tuple path(genome), path(fai)
    output: path 'star_index'
    script:
    """
    L=\$(awk '{s+=\$2} END{print s}' ${fai})
    nb=\$(python3 -c "import math; print(min(14, int(math.log2(\$L)/2 - 1)))")
    mkdir star_index
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} --genomeDir star_index \\
         --genomeFastaFiles ${genome} --genomeSAindexNbases \$nb --outFileNamePrefix star_index/
    """
}

process STAR_ALIGN {
    tag "${sample}"
    label 'medium'
    publishDir "${params.outdir}/star", mode: 'copy', pattern: '*.{bam,bai,Log.final.out,SJ.out.tab}'
    input:
    tuple val(sample), val(strand), path(r1), path(r2)
    path index
    output:
    tuple val(sample), val(strand), path("${sample}.bam"), path("${sample}.bam.bai"), emit: bam
    path "${sample}.Log.final.out", emit: log
    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir ${index} --readFilesIn ${r1} ${r2} \\
         --readFilesCommand ${r1.name.endsWith('.gz') ? 'zcat' : 'cat'} \\
         --outFileNamePrefix ${sample}. --outSAMtype BAM SortedByCoordinate \\
         --limitBAMsortRAM ${(task.memory.toBytes() * 0.5) as long} \\
         --outSAMattributes NH HI AS nM NM MD --outSAMattrRGline ID:${sample} SM:${sample} \\
         --outSAMunmapped None --twopassMode Basic --outFilterType BySJout \\
         --alignIntronMax ${params.intron_max} --alignMatesGapMax ${params.intron_max} \\
         --outTmpDir \$PWD/star_tmp
    mv ${sample}.Aligned.sortedByCoord.out.bam ${sample}.bam
    samtools index -@ ${task.cpus} ${sample}.bam
    rm -rf ${sample}._STARgenome ${sample}._STARpass1
    """
}

process BAM_PREP {
    // validate a user BAM (coordinate-sorted, paired) and index it
    tag "${sample}"
    label 'small'
    input:  tuple val(sample), val(strand), path(bam)
    output: tuple val(sample), val(strand), path(bam), path("${bam}.bai")
    script:
    """
    samtools view -H ${bam} | grep -q '^@HD.*SO:coordinate' || { echo "${sample}: BAM is not coordinate-sorted" >&2; exit 1; }
    np=\$(samtools view -c -f 1 ${bam} 2>/dev/null | head -1)
    [ "\$np" -gt 0 ] || { echo "${sample}: no paired reads; pipeline needs stranded paired-end data" >&2; exit 1; }
    samtools index ${bam}
    """
}

process JUNCTIONS {
    tag "${sample}"
    label 'small'
    input:  tuple val(sample), val(strand), path(bam), path(bai)
    output: path "${sample}.junctions.tsv"
    script:
    """
    samtools view -q ${params.min_mapq} -f 2 -F 0xF0C ${bam} \\
      | awk -F'\\t' '\$6 ~ /N/' \\
      | bam_junctions.py --strandedness ${strand} > ${sample}.junctions.tsv
    """
}

process STRAND_COVERAGE {
    // pooled coverage of one transcript strand for all libraries of one strandedness
    tag "${strand}:${tstrand}"
    label 'coverage'
    input:  tuple val(strand), val(tstrand), path(bams), path(bais)
    output: tuple val(tstrand), path("cov.${strand}.${tstrand}.wig")
    script:
    // reverse (dUTP): transcript strand = R2 strand = opposite of R1; forward: = R1 strand
    def plusRev = '(flag.read1 && flag.reverse) || (flag.read2 && !flag.reverse)'
    def plusFwd = '(flag.read1 && !flag.reverse) || (flag.read2 && flag.reverse)'
    def expr = (tstrand == 'plus') == (strand == 'reverse') ? plusRev : plusFwd
    def dflag = params.dedup_overlap ? '-s' : ''
    """
    out=cov.${strand}.${tstrand}.wig
    echo "track type=wiggle_0 name=${strand}_${tstrand}" > \$out
    # one seqid at a time: samtools depth -s buffers reads awaiting mates, which grows without bound
    # over a whole-genome stream (seen as an OOM at 16 GB on 42 yeast libraries)
    for c in \$(samtools view -H ${bams[0]} | awk -F'\\t' '\$1=="@SQ"{sub("SN:","",\$2); print \$2}'); do
      echo "variableStep chrom=\$c" >> \$out
      samtools merge -@ ${task.cpus} -u -R "\$c" -o - ${bams} \\
        | samtools view -u -q ${params.min_mapq} -f 2 -F 0xF0C -e '${expr}' - \\
        | samtools depth ${dflag} - \\
        | awk '{print \$2, \$3}' >> \$out
    done
    """
}

process SUM_COVERAGE {
    tag "${tstrand}"
    label 'small_mem'
    publishDir "${params.outdir}/utr/coverage", mode: 'copy'
    input:  tuple val(tstrand), path(wigs), path(fai)
    output: path "cov.${tstrand}.wig"
    script: "sum_wig.py ${fai} cov.${tstrand}.wig ${wigs}"
}

process PREP_INPUTS {
    label 'small_mem'
    publishDir "${params.outdir}/utr", mode: 'copy'
    input:
    tuple path(genome), path(fai)
    path gff
    path junctions
    output:
    path 'prep/*', emit: all
    path 'prep_report.txt', emit: report
    script:
    def ex = params.exclude_chroms ? "--exclude-chroms '${params.exclude_chroms}'" : ''
    """
    make_utr_inputs.py --gff ${gff} --genome ${genome} --junctions ${junctions} --outdir prep \\
        --min-intron-reads ${params.min_intron_reads} ${ex}
    cp prep/prep_report.txt .
    """
}

process UTRRNASEQ {
    tag "${label}:${tstrand}"
    label 'small_mem'
    publishDir "${params.outdir}/utr/utrrnaseq", mode: 'copy'
    input:
    tuple val(label), val(smooth), val(zero), val(tstrand), path(wig)
    tuple path(genome), path(fai)
    path prep
    output: tuple val(label), path("UTRs.${label}.${tstrand}.gff")
    script:
    """
    utrrnaseq --in-scaffold-file ${genome} --in-coding-region-file codons.${tstrand}.gtf \\
      --in-intron-file introns.gff --in-wiggle-file ${wig} --out-file-name UTRs.${label}.${tstrand}.gff \\
      --read-length ${params.read_len} --limit ${params.utr_limit} --smoothing-window-size ${smooth} \\
      --minimum-average-coverage 10 ${zero ? '--zero-coverage' : ''} > utrrnaseq.${label}.${tstrand}.log 2>&1
    """
}

process THRESHOLD_UTR {
    label 'calls'
    publishDir "${params.outdir}/utr", mode: 'copy'
    input:
    tuple path(genome), path(fai)
    path prep
    path wigs
    path ref
    path utrr
    val utrr_labels
    output:
    path 'UTRs.threshold.gff', emit: gff
    path 'utr_calls.tsv', emit: calls
    path 'calibration.tsv'
    path 'utr_summary.txt'
    script:
    def r = ref.name != 'NO_REF' ? "--ref ${ref} --ref-id-strip '${params.ref_id_strip}'" : ''
    def calibrate = ref.name != 'NO_REF' && !params.fix_f
    def f5 = calibrate ? '' : "--f5 ${params.f5}"
    def f3 = calibrate ? '' : "--f3 ${params.f3}"
    def u = utrr_labels.collect { l -> "--utrrnaseq ${l}=UTRs.${l}" }.join(' ')
    """
    threshold_utr.py --utrdir . --genome ${genome} ${r} ${f5} ${f3} ${u} \\
      --limit ${params.utr_limit} --smooth ${params.smooth} --floor ${params.floor} \\
      --min-cds-cov ${params.min_cds_cov} --min-junction-reads ${params.min_intron_reads} \\
      --junction-frac ${params.junction_frac} --seed ${params.seed} --summary utr_summary.txt
    """
}

process MERGE_GFF {
    label 'small_mem'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path gff
    path calls
    path utrgff
    output:
    path "${params.out_prefix}.utr.gff3", emit: gff
    path "${params.out_prefix}.utr.gff3.report.txt"
    path "${params.out_prefix}.utr.gffread_check.txt"
    script:
    """
    mkdir utr && cp ${calls} ${utrgff} utr/
    merge_utrs.py --gff ${gff} --utrdir utr --out ${params.out_prefix}.utr.gff3 ${params.keep_capped ? '--keep-capped' : ''}
    gffread -E ${params.out_prefix}.utr.gff3 -o /dev/null > ${params.out_prefix}.utr.gffread_check.txt 2>&1
    if ${params.strict_gff_check} && grep -qiE 'error|warning' ${params.out_prefix}.utr.gffread_check.txt; then
      echo "gffread reported problems:" >&2; cat ${params.out_prefix}.utr.gffread_check.txt >&2; exit 1
    fi
    """
}

process SALMON_INDEX {
    label 'medium'
    publishDir "${params.outdir}/salmon", mode: 'copy', pattern: '{transcripts.fa,tx2gene.tsv}'
    input:
    tuple path(genome), path(fai)
    path gff
    output:
    path 'index', emit: index
    path 'tx2gene.tsv', emit: tx2gene
    path 'transcripts.fa'
    script:
    """
    gffread -w transcripts.fa -g ${genome} ${gff}
    tx2gene.py ${gff} transcripts.fa > tx2gene.tsv
    cut -f1 ${fai} > decoys.txt
    cat transcripts.fa ${genome} > gentrome.fa
    mkdir -p idx_tmp
    salmon index -t gentrome.fa -d decoys.txt -i index -k ${params.salmon_k} -p ${task.cpus} --tmpdir idx_tmp
    rm -f gentrome.fa
    """
}

process SALMON_QUANT {
    tag "${sample}"
    label 'medium'
    publishDir "${params.outdir}/salmon/quant", mode: 'copy'
    input:
    tuple val(sample), val(strand), path(r1), path(r2)
    path index
    output: path "${sample}"
    script:
    """
    salmon quant -i ${index} -l ${params.salmon_libtype} -1 ${r1} -2 ${r2} \\
      ${params.salmon_bias ? '--gcBias --seqBias' : ''} -p ${task.cpus} -o ${sample} --radScratchDir \$PWD
    """
}

// ---------------------------------------------------------------- workflow
workflow {
    if (params.help) { helpMessage(); return }
    checkParams()
    def rows = readSamplesheet(params.samplesheet as String)
    def has_fq = rows.any { r -> r.fastq_1 }
    if (params.run_salmon && !has_fq) error "--run_salmon needs FASTQ rows"
    def ref_utr = refUtr()

    genome = FAIDX(channel.value(file(params.genome, checkIfExists: true)))
    gff    = channel.value(file(params.gff, checkIfExists: true))

    fq_rows  = channel.fromList(rows.findAll { r -> r.fastq_1 })
                      .map { r -> tuple(r.sample, r.strandedness, file(r.fastq_1, checkIfExists: true), file(r.fastq_2, checkIfExists: true)) }
    bam_rows = channel.fromList(rows.findAll { r -> r.bam })
                      .map { r -> tuple(r.sample, r.strandedness, file(r.bam, checkIfExists: true)) }

    def star_bams = channel.empty()
    if (has_fq) {
        index = params.star_index ? channel.value(file(params.star_index, checkIfExists: true)) : STAR_INDEX(genome)
        star_bams = STAR_ALIGN(fq_rows, index).bam
    }
    bams = star_bams.mix(BAM_PREP(bam_rows))

    junc = JUNCTIONS(bams).collect()

    cov_in = bams
        .map { _s, st, b, i -> tuple(st, b, i) }
        .groupTuple()
        .flatMap { st, bs, is -> ['plus', 'minus'].collect { ts -> tuple(st, ts, bs, is) } }
    wigs = STRAND_COVERAGE(cov_in)
        .groupTuple()
        .combine(genome.map { _g, f -> f })
    cov = SUM_COVERAGE(wigs)

    prep = PREP_INPUTS(genome, gff, junc).all

    // optional utrrnaseq comparison runs (tool default window, 50 bp window)
    def utrr_files  = channel.value([])
    def utrr_labels = []
    if (params.run_utrrnaseq) {
        runs = channel.fromList([['w150', 150, false], ['w50', 50, false]])
        cov_by_strand = cov.map { w -> tuple(w.name.contains('.plus.') ? 'plus' : 'minus', w) }
        u_in = runs.combine(cov_by_strand).map { l, sm, z, ts, w -> tuple(l, sm, z, ts, w) }
        utrr_files  = UTRRNASEQ(u_in, genome, prep).map { _l, f -> f }.collect()
        utrr_labels = ['w150', 'w50']
    }

    ref = ref_utr ? channel.value(file(ref_utr, checkIfExists: true)) : channel.value(file("${projectDir}/assets/NO_REF"))
    calls = THRESHOLD_UTR(genome, prep, cov.collect(), ref, utrr_files, channel.value(utrr_labels))

    merged = MERGE_GFF(gff, calls.calls, calls.gff).gff

    if (params.run_salmon) {
        sidx = SALMON_INDEX(genome, merged)
        SALMON_QUANT(fq_rows, sidx.index)
    }
}
