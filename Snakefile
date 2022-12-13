# 2022 Alex Kim
# Maintainer: Alex Kim
# Email: alex.kim@agresearch.co.nz


configfile: 

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
SAMPLES, = glob_wildcards('fastq/{samples}.fastq.gz')

# sanity check
print("Found: ")
for WLDCRD in SAMPLES:
    print(WLDCRD)
print("")

rule all:
    input: 
        # kraken2 report 
        # multiqc report 
        'results/ReadsMultiQCReport.html'
        

rule fastqc:
    # quality control 
    input:
        fastq = expand('fastq/{samples}.fastq.gz', samples = SAMPLES)
    output: 
        html = 'results/fastqc/{samples}_fastqc.html',
        zip = 'results/fastqc/{samples}_fastqc.zip'
    conda: 
        'envs/fastqc.yaml'
    threads: 1 
    message: 
        'Running QC on reads: {wildcards.samples}\n'
    shell:
        '-fastqc '
        '-o results/fastqc/ '
        '-q '
        '-t {threads} '
        '{input.fastq}'


rule multiqc:
    # reporting tool 

rule kneaddata: 
    # quality control 

rule kraken2:
    # taxonomic profiling 

rule braken:
    # compute abundance 

#  merge kraken2 reports using utility scripts 


