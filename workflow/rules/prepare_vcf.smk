# Create population-specific sample lists. Each list contains target samples
# followed by the common outgroup required by SPrime.
rule samples:
    input: qc='qc/validation.tsv', samples=C['samples']
    output: samples='work/{population}/samples.txt', outgroup='work/{population}/outgroup.txt'
    params: stage='samples', code=CODE
    log: 'logs/{population}/samples.log'
    script: TASK

# Filter one population/chromosome to biallelic SNPs and retain FORMAT/GT only.
rule subset:
    input: samples='work/{population}/samples.txt', vcf=lambda w: C['vcf'].format(chrom=w.chrom)
    output: vcf='work/{population}/chr{chrom}.vcf.gz', index='work/{population}/chr{chrom}.vcf.gz.csi'
    params: stage='subset', code=CODE
    resources: mem_mb=2000
    log: 'logs/{population}/subset.{chrom}.log'
    script: TASK

# SPrime receives a single ordered autosomal VCF for each population.
rule concat:
    input: vcfs=lambda w: expand('work/{population}/chr{chrom}.vcf.gz', population=w.population, chrom=CH)
    output: vcf='work/{population}/all.auto.vcf.gz'
    params: stage='concat', code=CODE
    log: 'logs/{population}/concat.log'
    script: TASK
