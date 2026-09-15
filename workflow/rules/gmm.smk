# Recount match rates after the strict SPrime SCORE filter used by pre_data.r.
# No physical segment-length filter is applied.
rule prepare_gmm_base:
    input: lambda w: expand('archaic_match/{population}/{ref}/chr{chrom}.mscore', population=w.population, ref=list(R), chrom=CH)
    output: 'work/{population}/gmm_base.tsv'
    params: stage='gmm_base', score_gt=C['gmm']['score_gt'], code=CODE
    resources: mem_mb=8000
    log: 'logs/{population}/gmm_base.log'
    script: TASK

# Apply callable-allele and Nean/Deni match-rate eligibility conditions.
rule gmm_input:
    input: 'work/{population}/gmm_base.tsv'
    output: 'gmm/input/{population}.tsv'
    params: stage='gmm_input', settings=C['gmm'], code=CODE
    log: 'logs/{population}/gmm_input.log'
    script: TASK

# Fit 1/2/3 components independently for each population and target Deni ref.
rule fit_gmm:
    input: 'gmm/input/{population}.tsv'
    output: selection='gmm/{population}/{ref}.selection.tsv', components='gmm/{population}/{ref}.components.tsv', plot='gmm/plots/{population}.{ref}.png'
    params: stage='gmm', settings=C['gmm'], code=CODE, populations=P
    resources: mem_mb=2000
    log: 'logs/{population}/gmm.{ref}.log'
    script: TASK

# Produce project-level model-selection and component parameter tables.
rule collect_gmm:
    input:
        selection=expand('gmm/{population}/{ref}.selection.tsv', population=P, ref=C['gmm']['target_references']),
        components=expand('gmm/{population}/{ref}.components.tsv', population=P, ref=C['gmm']['target_references']),
        plots=expand('gmm/plots/{population}.{ref}.png', population=P, ref=C['gmm']['target_references'])
    output: selection='gmm/model_selection.tsv', components='gmm/components.tsv'
    params: stage='collect_gmm', code=CODE
    log: 'logs/collect_gmm.log'
    script: TASK
