# Call carried SPrime-marker runs for both haplotypes of every target sample.
rule individual_chr:
    input: vcf='work/{population}/chr{chrom}.vcf.gz', index='work/{population}/chr{chrom}.vcf.gz.csi', score='sprime/{population}/chr{chrom}.score', classes='classification/{population}.tsv', samples=C['samples']
    output: 'work/{population}/individual.{chrom}.tsv'
    params: stage='individual', settings=C['individual'], code=SIG['individual']
    resources: mem_mb=C['resources']['individual_mem_mb']
    log: 'logs/{population}/individual.{chrom}.log'
    script: TASK

# Merge chromosomes and provide convenient class-specific population files.
rule individual_pop:
    input: lambda w: expand('work/{population}/individual.{chrom}.tsv', population=w.population, chrom=CH)
    output:
        all='individual_calls/by_population/{population}/all.tsv.gz',
        nean='individual_calls/by_population/{population}/neanderthal.tsv.gz',
        deni='individual_calls/by_population/{population}/denisovan.tsv.gz',
        ambiguous='individual_calls/by_population/{population}/ambiguous.tsv.gz'
    params: stage='individual_pop', code=SIG['individual_pop']
    log: 'logs/{population}/individual_collect.log'
    script: TASK

# Final all-population individual-call table.
rule collect_individual:
    input: expand('individual_calls/by_population/{population}/all.tsv.gz', population=P)
    output: 'individual_calls/all_individual_calls.tsv.gz'
    params: stage='collect', code=SIG['collect']
    log: 'logs/individual_collect.log'
    script: TASK

# Summarize haplotype-specific calls per target individual and class.
rule individual_summary:
    input:
        calls=expand('individual_calls/by_population/{population}/all.tsv.gz', population=P),
        samples=C['samples']
    output: 'individual_calls/individual_summary.tsv'
    params: stage='individual_summary', code=SIG['individual_summary']
    log: 'logs/individual_summary.log'
    script: TASK
