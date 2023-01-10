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
        expand('results/kraken2GTDB/{samples}.GTDB.k2report', samples = SAMPLES),
        # bracken outputs
        # expand('results/bracken/{samples}.bracken', samples = SAMPLES),
        # multiqc reports (raw data and knead data)
        'results/ReadsMultiQCReportRawData.html',
        'results/ReadsMultiQCReportKneadData.html',
        expand('results/humann3/{samples}_genefamilies.tsv', samples = SAMPLES),
        expand('results/humann3protein/{samples}_genefamilies.tsv', samples = SAMPLES)
        

rule merge:
    input: 
        read1 = 'fastq/{samples}_R1_001.fastq.gz',
        read2 = 'fastq/{samples}_R2_001.fastq.gz'
    output:
        mergedReads = 'results/mergedReads/{samples}.fastq.gz'
    shell:
        'cat {input.read1} {input.read2} > {output.mergedReads}'


rule fastqc:
    # quality control 
    input:
        fastq = 'results/mergedReads/{samples}.fastq.gz'
    output: 
        html = 'results/fastqc/{samples}_fastqc.html',
        zip = 'results/fastqc/{samples}_fastqc.zip'
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
        fastqc = expand('results/fastqc/{samples}_fastqc.zip', samples = SAMPLES)
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
        fastq = 'results/mergedReads/{samples}.fastq.gz'
    output: 
        # trim adapters ...
        trimReads = temp('results/kneaddata/{samples}_kneaddata.trimmed.fastq'),
        # trim repetitive sequences 
        trfReads = temp('results/kneaddata/{samples}_kneaddata.repeats.removed.fastq'),
        # trim host DNA 
        ovineReads = temp('results/kneaddata/{samples}_kneaddata_ARS_UI_Ramb_v2_bowtie2_contam.fastq'),
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
        '-db /bifo/scratch//2022-AK-MBIE-Rumen-MG/ref/ARS_UI_Ramb_v2 '
        '-db /bifo/scratch//2022-AK-MBIE-Rumen-MG/ref/SILVA_128_LSUParc_SSUParc_ribosomal_RNA '
        '-o results/kneaddata && '
        'seqkit stats -j 12 -a results/kneaddata/{wildcards.samples}*.fastq > {output.readStats}'


rule fastqcKDR: 
    input: 
        fastqc = 'results/kneaddata/{samples}_kneaddata.fastq'
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


rule kraken2GTDB:
    # taxonomic profiling 
    input:
        KDRs = 'results/kneaddata/{samples}_kneaddata.fastq'
    output: 
        k2OutputGTDB = 'results/kraken2GTDB/{samples}.GTDB.kraken2',
        k2ReportGTDB = 'results/kraken2GTDB/{samples}.GTDB.k2report'
    log:
        'logs/kraken2GTDB/{samples}.kraken2.GTDB.log'
    conda:
        'envs/kraken2.yaml'
    threads: 20
    resources: 
        mem_gb=360,
        partition="inv-bigmem-fast"
    shell:
        'kraken2 '
        '--use-names '
        '--db /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB '
        '-t {threads} '
        '--report {output.k2ReportGTDB} '
        '--report-minimizer-data '
        '{input.KDRs} > {output.k2OutputGTDB}'


rule bracken:
    # compute abundance 
    input:
        k2ReportGTDB = 'results/kraken2GTDB/{samples}.GTDB.k2report'
    output:
        'results/bracken/{samples}.bracken',
        'results/bracken/{samples}.breport'
    log:
        'logs/bracken/{samples}.bracken.log'
    conda:
        'envs/bracken.yaml'
    threads: 8
    shell: 
        'bracken '
        '--log {log} '
        '-d /dataset/2022-BJP-GTDB/scratch/2022-BJP-GTDB/kraken/GTDB '
        '-i results/kraken2GTDB/{wildcards.samples}.GTDB.k2report '
        '-o results/bracken/{wildcards.samples}.bracken '
        '-w results/bracken/{wildcards.samples}.breport '
        '-r 240 ' # average read length
        '-l S '  # species
        '-t 10' # remove low abundance species (noise)  


rule humann3:
    # functional profiling
    input:
        KDRs = rules.kneaddata.output.clnReads
    output:
        geneFamilies = 'results/humann3/{samples}_genefamilies.tsv',
        pathwaysCoverage = 'results/humann3/{samples}_pathcoverage.tsv',
        pathways = 'results/humann3/{samples}_pathabundance.tsv'
    log:
        'logs/humann3/{samples}.humann3.log'
    conda:
        'envs/humann3.yaml'
    threads: 10
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann3 ' 
        '--threads {threads} '
        '--bypass-nucleotide-index '
        '--search-mode uniref50 '
        '--nucleotide-database /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/humann3Struo2/uniref50 '
        '--protein-database /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/humann3Struo2/uniref50/protein_database '
        '--input-format fastq '
        '--output results/humann3 '
        '--input {input.KDRs} '
        '--output-basename {samples} '
        '--o-log {log} '

rule humann3protein:
    # functional profiling
    input:
        KDRs = rules.kneaddata.output.clnReads
    output:
        geneFamilies = 'results/humann3protein/{samples}_genefamilies.tsv',
        pathwaysCoverage = 'results/humann3protein/{samples}_pathcoverage.tsv',
        pathways = 'results/humann3protein/{samples}_pathabundance.tsv'
    log:
        'logs/humann3protein/{samples}.humann3.log'
    conda:
        'envs/humann3.yaml'
    threads: 10
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann3 ' 
        '--threads {threads} '
        '--bypass-nucleotide-search '
        '--search-mode uniref50 '
        '--protein-database /bifo/scratch/2022-BJP-GTDB/2022-BJP-GTDB/humann3Struo2/uniref50/protein_database '
        '--input-format fastq '
        '--output results/humann3protein '
        '--input {input.KDRs} '
        '--output-basename {samples} '
        '--o-log {log} '
