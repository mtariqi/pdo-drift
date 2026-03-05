# ─────────────────────────────────────────────────────────────
#  WES Pipeline — Somatic Variant Calling
#  Tools: fastp → bwa-mem2 → GATK4 MuTect2 → CNVKit → PyClone-VI
# ─────────────────────────────────────────────────────────────
import pandas as pd
from pathlib import Path

configfile: "../../config/config.yaml"

SAMPLES = pd.read_csv(config["samples"], sep="\t")
TUMOR_IDS = SAMPLES[SAMPLES["condition"] == "fresh"]["sample_id"].tolist()

rule all:
    input:
        expand("results/variants/{sample}.filtered.vcf.gz", sample=TUMOR_IDS),
        expand("results/cnv/{sample}.cnr", sample=TUMOR_IDS),
        expand("results/clonal/{sample}.pyclone.tsv", sample=TUMOR_IDS),
        "results/qc/multiqc_wes_report.html"

# ── Step 1: Read trimming & QC ────────────────────────────────
rule fastp_trim:
    input:
        r1="data/raw/{sample}_R1.fastq.gz",
        r2="data/raw/{sample}_R2.fastq.gz"
    output:
        r1="data/trimmed/{sample}_R1.fastq.gz",
        r2="data/trimmed/{sample}_R2.fastq.gz",
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json"
    threads: config["resources"]["default_threads"]
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --html {output.html} --json {output.json} \
            --detect_adapter_for_pe \
            --correction \
            --thread {threads} \
            --qualified_quality_phred 20 \
            --length_required 50
        """

# ── Step 2: Alignment with bwa-mem2 ──────────────────────────
rule bwa_align:
    input:
        r1="data/trimmed/{sample}_R1.fastq.gz",
        r2="data/trimmed/{sample}_R2.fastq.gz",
        idx=config["reference"]["genome_index_bwa"]
    output:
        bam=temp("data/aligned/{sample}.raw.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA\tLB:lib1"
    threads: config["resources"]["alignment_threads"]
    shell:
        """
        bwa-mem2 mem \
            -t {threads} \
            -R "{params.rg}" \
            {input.idx}/hg38.fa \
            {input.r1} {input.r2} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# ── Step 3: Mark duplicates ───────────────────────────────────
rule mark_duplicates:
    input:
        "data/aligned/{sample}.raw.bam"
    output:
        bam="data/aligned/{sample}.markdup.bam",
        metrics="results/qc/picard/{sample}.dup_metrics.txt"
    shell:
        """
        picard MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.metrics} \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT
        """

# ── Step 4: Base Quality Score Recalibration ─────────────────
rule bqsr:
    input:
        bam="data/aligned/{sample}.markdup.bam",
        ref=config["reference"]["genome"],
        dbsnp=config["reference"]["dbsnp"]
    output:
        bam="data/aligned/{sample}.bqsr.bam",
        table="data/aligned/{sample}.recal.table"
    params:
        mem=config["resources"]["gatk_mem"]
    shell:
        """
        gatk --java-options "-Xmx{params.mem}" BaseRecalibrator \
            -I {input.bam} -R {input.ref} \
            --known-sites {input.dbsnp} \
            -O {output.table}

        gatk --java-options "-Xmx{params.mem}" ApplyBQSR \
            -I {input.bam} -R {input.ref} \
            --bqsr-recal-file {output.table} \
            -O {output.bam}
        """

# ── Step 5: Somatic variant calling (MuTect2) ────────────────
rule mutect2:
    input:
        tumor="data/aligned/{sample}.bqsr.bam",
        normal=lambda wc: "data/aligned/{}_normal.bqsr.bam".format(
            SAMPLES.loc[SAMPLES["sample_id"]==wc.sample, "patient_id"].values[0]
        ),
        ref=config["reference"]["genome"],
        gnomad=config["reference"]["gnomad"]
    output:
        vcf="results/variants/{sample}.raw.vcf.gz",
        f1r2="results/variants/{sample}.f1r2.tar.gz",
        stats="results/variants/{sample}.vcf.stats"
    params:
        mem=config["resources"]["gatk_mem"]
    shell:
        """
        gatk --java-options "-Xmx{params.mem}" Mutect2 \
            -R {input.ref} \
            -I {input.tumor} \
            -I {input.normal} \
            --normal-sample $(basename {input.normal} .bqsr.bam) \
            --germline-resource {input.gnomad} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf}
        """

# ── Step 6: Filter variants ───────────────────────────────────
rule filter_variants:
    input:
        vcf="results/variants/{sample}.raw.vcf.gz",
        stats="results/variants/{sample}.vcf.stats",
        ref=config["reference"]["genome"]
    output:
        vcf="results/variants/{sample}.filtered.vcf.gz"
    params:
        mem=config["resources"]["gatk_mem"],
        min_af=config["wes"]["min_allele_freq"]
    shell:
        """
        gatk --java-options "-Xmx{params.mem}" FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            --stats {input.stats} \
            --min-allele-fraction {params.min_af} \
            -O {output.vcf}
        """

# ── Step 7: CNV calling (CNVKit) ──────────────────────────────
rule cnvkit:
    input:
        tumor="data/aligned/{sample}.bqsr.bam",
        targets=config["reference"]["exome_targets"],
        ref=config["reference"]["genome"]
    output:
        cnr="results/cnv/{sample}.cnr",
        cns="results/cnv/{sample}.cns"
    threads: config["resources"]["default_threads"]
    shell:
        """
        cnvkit.py batch {input.tumor} \
            --method hybrid \
            --targets {input.targets} \
            --fasta {input.ref} \
            --output-dir results/cnv/ \
            -p {threads}
        """

# ── Step 8: Clonal evolution (PyClone-VI) ────────────────────
rule pyclone:
    input:
        vcf="results/variants/{sample}.filtered.vcf.gz",
        cnv="results/cnv/{sample}.cns"
    output:
        "results/clonal/{sample}.pyclone.tsv"
    params:
        num_clusters=10,
        num_restarts=20
    shell:
        """
        python ../../src/genomics/run_pyclone.py \
            --vcf {input.vcf} \
            --cnv {input.cnv} \
            --output {output} \
            --num-clusters {params.num_clusters} \
            --num-restarts {params.num_restarts}
        """

# ── Step 9: Aggregate QC ──────────────────────────────────────
rule multiqc:
    input:
        expand("results/qc/fastp/{sample}.json", sample=TUMOR_IDS),
        expand("results/qc/picard/{sample}.dup_metrics.txt", sample=TUMOR_IDS)
    output:
        "results/qc/multiqc_wes_report.html"
    shell:
        "multiqc results/qc/ -o results/qc/ -n multiqc_wes_report"
