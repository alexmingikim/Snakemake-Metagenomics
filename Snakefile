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

READS = ['_R1_001', '_R2_001']

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
        

rule fastqcR1:
    # quality control 
    input:
        # QC for files individually
        fastq = 'fastq/{samples}_R1_001.fastq.gz'
    output: 
        html = 'results/fastqc/{samples}_R1_001_fastqc.html',
        zip = 'results/fastqc/{samples}_R1_001_fastqc.zip'
    conda: 
        'envs/fastqc.yaml'
    threads: 1 
    message: 
        'Running quality checks on reads: {wildcards.samples}_R1_001\n'
    shell:
        'fastqc '
        '-o results/fastqc/ '
        '-q ' # suppress progress messages; only report errors 
        '-t {threads} '
        '{input.fastq}'


rule fastqcR2:
    # quality control 
    input:
        # QC for files individually
        fastq = 'fastq/{samples}_R2_001.fastq.gz'
    output: 
        html = 'results/fastqc/{samples}_R2_001_fastqc.html',
        zip = 'results/fastqc/{samples}_R2_001_fastqc.zip'
    conda: 
        'envs/fastqc.yaml'
    threads: 1 
    message: 
        'Running quality checks on reads: {wildcards.samples}_R2_001\n'
    shell:
        'fastqc '
        '-o results/fastqc/ '
        '-q ' # suppress progress messages; only report errors 
        '-t {threads} '
        '{input.fastq}'

rule multiqc:
    # reporting tool 
    input: 
        fastqc = expand('results/fastqc/{samples}{reads}_fastqc.zip', samples = SAMPLES, reads = READS)
    output: 
        multiqc = 'results/ReadsMultiQCReportRawData.html'
    conda: 
        "envs/multiqc.yaml"
    shell:
        'multiqc '
        '-n results/ReadsMultiQCReportRawData '
        '-s ' # to not clean sample names 
        '-f ' # overwrite existing reports 
        '--interactive ' # interactive plots 
        '{input.fastqc}'


rule kneaddata: 
    # quality control - separate bacterial reads from contaminant reads (host, bacterial 16S sequences) 
    input: 
        ### two inputs for paired end reads? ###
        fastq1 = 'fastq/{samples}_R1_001.fastq.gz', 
        fastq2 = 'fastq/{samples}_R2_001.fastq.gz'
    output: 
        """
        # trim adapters ...
        trimReads = temp('results/kneaddata/{samples}_kneaddata.trimmed.fastq'),
        # trim repetitive sequences 
        trfReads = temp('results/kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        # trim host DNA 
        ovineReads = temp('results/kneaddata/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam.fastq'),
        # trim 16S rRNA
        silvaReads = temp('results/kneaddata/{samples}_kneaddata_SILVA_128_LSUParc_SSUParc_ribosomal_RNA_bowtie2_contam.fastq'),
        # filtered reads 
        clnReads1 = 'results/kneaddata/{samples}_R1_001_kneaddata_paired_1.fastq'
        clnReads2 = 'results/kneaddata/{samples}_R1_
        # summary?
        # readStats = 'results/kneaddata/{samples_short}.read.stats.txt'
        """
        # reads from R1 identified as NOT belonging to any reference databases 
        clnReadsR1 = 'results/kneaddata/{sample}_R1_001_kneaddata_paired_1.fastq'
        # reads from R2 identified as NOT belonging to any reference databases 
        clnReadsR2 = 'results/kneaddata/{sample}_R2_001_kneaddata_paired_2.fastq'

        # cases when one of the reads do not pass quality filtering  
        unmatchedR1 = temp('results/kneaddata/{sample}_R1_001_kneaddata_unmatched_1.fastq')
        unmatchedR2 = temp('results/kneaddata/{sample}_R2_001_kneaddata_unmatched_2.fastq')
        
    conda: 
        'envs/biobakery.yaml'
    log:
        'logs/kneaddata/{samples}.kneaddata.log'
    threads: 4
    message:
        'kneaddata: {wildcards.samples}\n'
    shell:
        'kneaddata '
        '--input1 {input.fastq1} '
        '--input2 {input.fastq2} '
        '-t {threads} '
        '--log-level INFO '
        '--log {log} '
        '--trimmomatic /home/kima/conda-envs/biobakery/share/trimmomatic-0.39-2 ' 
        '--sequencer-source TruSeq3 ' # to identify correct adapter sequences
        '-db ref/ARS_UI_Ramb_v2 '
        '-db ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA '
        '-o results/kneaddata'
        # 'seqkit stats -j 12 -a results/kneaddata/{wildcards.samples_long}*.fastq > {output.readStats}'


rule fastqcKDR1: 
    input: 
        fastqc = rules.kneaddata.output.clnReadsR1
    output:
        'results/fastqcKDR/{samples}_R1_001_kneaddata_fastqc.zip'
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

    
rule fastqcKDR2: 
    input: 
        fastqc = rules.kneaddata.output.clnReadsR2
    output:
        'results/fastqcKDR/{samples}_R2_001_kneaddata_fastqc.zip'
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
        fastqc = expand('results/fastqcKDR/{samples}{reads}_kneaddata_fastqc.zip', samples = SAMPLES, reads = READS)
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
        KDRs = rules.kneaddata.output.cleanReads
    output: 
        k2Output = 'results/kraken2/{samples_long}.out.k2',
        k2Report = 'results/kraken2/{samples_long}.report.k2'
    log:
        'logs/{samples_long}.kraken2.log'
    conda:
        'envs/kraken2.yaml'
    threads: 12
    resources: mem_mb=180000
    shell:
        'kraken2 '
        '--db ref/kraken2 '
        '--report {output.k2Report} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2Output)'
"""

rule braken:
    # compute abundance 

#  merge kraken2 reports using utility scripts 