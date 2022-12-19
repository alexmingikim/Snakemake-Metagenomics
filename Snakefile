# 2022 Alex Kim
# Maintainer: Alex Kim
# Email: alex.kim@agresearch.co.nz


# configfile: 

import os

onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')
    print(f"Env TMPDIR={os.environ.get('TMPDIR', '<n/a>')}")

# define samples from data directory using wildcards
SAMPLES, = glob_wildcards('fastq/{samples}_R1_001.fastq.gz') 


# sanity check
print("Found: ")
for WLDCRD in SAMPLES:
    print(WLDCRD)
print("")

rule all:
    input: 
        # kraken2 report 
        # expand('results/kraken2/{samples}.report.k2', samples = SAMPLES), 
        # multiqc reports (raw data and knead data)
        'results/ReadsMultiQCReportRawData.html',
        'results/ReadsMultiQCReportKneadData.html'
        

rule merge:
    input: 
        read1 = 'fastq/{samples}_R1_001.fastq.gz',
        read2 = 'fastq/{samples}_R2_001.fastq.gz'
    output:
        mergedReads = 'fastq/mergedReads/{samples}_merged_fastq.gz'
    shell:
        'cat fastq/{wildcards.samples}_R1_001.fastq.gz fastq/{wildcards.samples}_R2_001.fastq.gz > fastq/mergedReads/{wildcards.samples}_merged_fastq.gz'


rule fastqc:
    # quality control 
    input:
        fastq = 'fastq/mergedReads/{samples}_merged_fastq.gz'
    output: 
        html = 'results/fastqc/{samples}_merged_fastqc.html',
        zip = 'results/fastqc/{samples}_merged_fastqc.zip'
    conda: 
        'envs/fastqc.yaml'
    threads: 1 
    message: 
        'Running quality checks on reads: {wildcards.samples}\n'
    shell:
        'fastqc '
        '-o results/fastqc/ '
        '-q ' 
        '-t {threads} '
        '{input.fastq}'


rule multiqc:
    # reporting tool 
    input: 
        fastqc = expand('results/fastqc/{samples}_merged_fastqc.zip', samples = SAMPLES)
    output: 
        multiqc = 'results/ReadsMultiQCReportRawData.html'
    conda: 
        'envs/multiqc.yaml'
    shell:
        'multiqc '
        '-n results/ReadsMultiQCReportRawData '
        '-s ' 
        '-f ' 
        '--interactive ' 
        '{input.fastqc}'


rule kneaddata: 
    # quality control - separate bacterial reads from contaminant reads (host, bacterial 16S sequences) 
    input: 
        fastq = rules.merge.output.mergedReads
    output: 
        # trim adapters ...
        trimReads = temp('results/kneaddata/{samples}_kneaddata.trimmed.fastq'),
        # trim repetitive sequences 
        trfReads = temp('results/kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        # trim host DNA 
        ovineReads = temp('results/kneaddata/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie_contam.fastq'),
        # trim 16S rRNA
        silvaReads = temp('results/kneaddata/{samples}_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq'),
        # filtered reads 
        clnReads = 'results/kneaddata/{samples}_kneaddata.fastq',
        # seqkit stats
        readStats = 'results/seqkit/{samples}.read.stats.txt'
    conda: 
        'envs/biobakery.yaml'
    log:
        'logs/kneaddata/{samples}.kneaddata.log'
    threads: 4
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--unpaired {input.fastq} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/kima/conda-envs/biobakery/share/trimmomatic-0.39-2 ' 
        '--sequencer-source TruSeq3 ' # to identify correct adapter sequences
        '-db ref/ARS_UI_Ramb_v2 '
        '-db ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA '
        '-o results/kneaddata && '
        'seqkit stats -j 12 -a results/kneaddata/{wildcards.samples}*.fastq > {output.readStats}'


rule fastqcKDR: 
    input: 
        fastqc = rules.kneaddata.output.clnReads
    output:
        zip = 'results/fastqcKDR/{samples}_kneaddata_fastqc.zip'
    conda: 
        'envs/fastqc.yaml'
    threads: 1
    message: 
        'Running quality checks on reads: {wildcards.samples}\n'
    shell: 
        'fastqc '
        '-o results/fastqcKDR/ '
        '-q '
        '-t {threads} '
        '{input.fastqc}'


rule multiQCKDRs: 
    input: 
        fastqc = expand('results/fastqcKDR/{samples}_kneaddata_fastqc.zip', samples = SAMPLES)
    output: 
        'results/ReadsMultiQCReportKneadData.html'
    conda: 
        'envs/multiqc.yaml'
    shell:
        'multiqc '
        '-n results/ReadsMultiQCReportKneadData '
        '-s '
        '-f '
        '--interactive '
        '{input.fastqc}'


"""
rule kraken2:
    # taxonomic profiling 
    input:
        reads1 = rules.kneaddata.output.clnReadsR1,
        reads2 = rules.kneaddata.output.clnReadsR2
    output: 
        k2Output = 'results/kraken2/{samples}.out.k2',
        k2Report = 'results/kraken2/{samples}.report.k2'
    log:
        'logs/{samples}.kraken2.log'
    conda:
        'envs/kraken2.yaml'
    threads: 12
    resources: 
        mem_gb=146,
        partition="inv-bigmem,inv-bigmem-fast"
    shell:
        'kraken2 '
        '--use-names '
        '--db ref/kraken2 '
        '-t {threads} '
        '--report {output.k2Report} '
        '--report-minimizer-data '
        '{input.reads1}'
"""

rule braken:
    # compute abundance 

#  merge kraken2 reports using utility scripts 