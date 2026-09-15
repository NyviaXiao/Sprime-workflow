# Detect candidate introgressed markers for one population/chromosome.
rule sprime:
    input: vcf='work/{population}/all.auto.vcf.gz', outgroup='work/{population}/outgroup.txt', jar=C['sprime_jar'], map=C['genetic_map']
    output: score='sprime/{population}/chr{chrom}.score'
    params: stage='sprime', minscore=C['sprime_minscore'], code=CODE
    resources: mem_mb=C['resources']['sprime_mem_mb']
    log: 'logs/{population}/sprime.{chrom}.log'
    script: TASK
