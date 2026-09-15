# Fail before expensive computation if sample/reference assumptions are broken.
rule validate:
    input:
        samples=C['samples'], genetic_map=C['genetic_map'], jar=C['sprime_jar'], binary=C['map_arch'],
        vcfs=sorted({C['vcf'].format(chrom=c) for c in CH}),
        archaic=sorted({r[k].format(chrom=c) for r in R.values() for c in CH for k in ('vcf','mask')})
    output: 'qc/validation.tsv', 'run_info.json'
    params: stage='validate', code=SIG['validate'], genome_build=C['genome_build'], chromosomes=CH
    log: 'logs/validate.log'
    script: TASK
