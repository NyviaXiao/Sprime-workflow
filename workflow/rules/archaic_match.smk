# Annotate one SPrime score file against one ancient individual. The reference
# dimension is fully driven by ``archaic_references`` in config.yaml.
rule archaic_match:
    input:
        score='sprime/{population}/chr{chrom}.score', binary=C['map_arch'],
        vcf=lambda w: R[w.ref]['vcf'].format(chrom=w.chrom),
        mask=lambda w: R[w.ref]['mask'].format(chrom=w.chrom)
    output: score='archaic_match/{population}/{ref}/chr{chrom}.mscore'
    params: stage='match', tag=lambda w: R[w.ref]['tag'], code=SIG['match']
    resources: mem_mb=C['resources']['map_arch_mem_mb']
    log: 'logs/{population}/match.{ref}.{chrom}.log'
    script: TASK
