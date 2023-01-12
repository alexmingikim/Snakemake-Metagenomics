# 2023 Alex Kim
# Maintainer: Alex Kim
# Email: alex.kim@agresearch.co.nz


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
SAMPLES, = glob_wildcards('KDRs/{samples}_kneaddata.fastq') 


# sanity check
print("Found: ")
for WLDCRD in SAMPLES:
    print(WLDCRD)
print("")


rule all:
    input: 
        expand('results/humann3Uniref50EC/{samples}_genefamilies.tsv', samples = SAMPLES),
        expand('results/humann3Uniref50Full/{samples}_genefamilies.tsv', samples = SAMPLES),
        expand('results/humann3Uniref90EC/{samples}_genefamilies.tsv', samples = SAMPLES)


rule humann3Uniref50EC:
    input:
        KDRs = 'KDRs/{samples}_kneaddata.fastq'
    output:
        geneFamilies = 'results/humann3Uniref50EC/{samples}_genefamilies.tsv',
        pathwaysCoverage = 'results/humann3Uniref50EC/{samples}_pathcoverage.tsv',
        pathways = 'results/humann3Uniref50EC/{samples}_pathabundance.tsv'
    log:
        'logs/humann3Uniref50EC/{samples}.humann3.log'
    conda:
        'envs/humann3.yaml'
    threads: 10
    resources: 
        mem_gb= lambda wildcards, attempt: attempt * 12,
        partition="inv-iranui",
        time="96:00:00"
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann3 ' 
        '--memory-use minimum '
        '--threads {threads} '
        '--bypass-nucleotide-search '
        '--search-mode uniref50 '
        '--protein-database /bifo/scratch/2022-AK-MBIE-Rumen-MG/ref/humann3/unirefECFilt '
        '--input-format fastq '
        '--output results/humann3Uniref50EC '
        '--input {input.KDRs} '
        '--output-basename {samples} '
        '--o-log {log}'

rule humann3Uniref50Full:
    input:
        KDRs = 'KDRs/{samples}_kneaddata.fastq'
    output:
        geneFamilies = 'results/humann3Uniref50Full/{samples}_genefamilies.tsv',
        pathwaysCoverage = 'results/humann3Uniref50Full/{samples}_pathcoverage.tsv',
        pathways = 'results/humann3Uniref50Full/{samples}_pathabundance.tsv'
    log:
        'logs/humann3Uniref50Full/{samples}.humann3.log'
    conda:
        'envs/humann3.yaml'
    threads: 10
    resources: 
        mem_gb= lambda wildcards, attempt: attempt * 12,
        partition="inv-iranui",
        time="96:00:00"
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann3 ' 
        '--memory-use minimum '
        '--threads {threads} '
        '--bypass-nucleotide-search '
        '--search-mode uniref50 '
        '--protein-database /bifo/scratch/2022-AK-MBIE-Rumen-MG/ref/humann3/uniref '
        '--input-format fastq '
        '--output results/humann3Uniref50Full '
        '--input {input.KDRs} '
        '--output-basename {samples} '
        '--o-log {log}'


rule humann3Uniref90EC:
    input:
        KDRs = 'KDRs/{samples}_kneaddata.fastq'
    output:
        geneFamilies = 'results/humann3Uniref90EC/{samples}_genefamilies.tsv',
        pathwaysCoverage = 'results/humann3Uniref90EC/{samples}_pathcoverage.tsv',
        pathways = 'results/humann3Uniref90EC/{samples}_pathabundance.tsv'
    log:
        'logs/humann3Uniref90EC/{samples}.humann3.log'
    conda:
        'envs/humann3.yaml'
    threads: 10
    resources: 
        mem_gb= lambda wildcards, attempt: attempt * 12,
        partition="inv-iranui",
        time="96:00:00"
    message:
        'humann3 profiling: {wildcards.samples}\n'
    shell:
        'humann3 ' 
        '--memory-use minimum '
        '--threads {threads} '
        '--bypass-nucleotide-search '
        '--search-mode uniref90 '
        '--protein-database /bifo/scratch/2022-AK-MBIE-Rumen-MG/ref/humann3/unirefECFilt '
        '--input-format fastq '
        '--output results/humann3Uniref90EC '
        '--input {input.KDRs} '
        '--output-basename {samples} '
        '--o-log {log}'