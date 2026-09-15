# Collapse marker-level mscore files into one row per SPrime segment, with a
# set of match/count columns for every configured ancient reference.
rule summarize:
    input: lambda w: expand('archaic_match/{population}/{ref}/chr{chrom}.mscore', population=w.population, ref=list(R), chrom=CH)
    output: summary='tables/{population}.wide.tsv'
    params: stage='summarize', summary=C['summary_min_callable'], code=CODE
    resources: mem_mb=8000
    log: 'logs/{population}/summary.log'
    script: TASK

# Combine population tables and emit both wide and tidy long representations.
rule collect_tables:
    input: expand('tables/{population}.wide.tsv', population=P)
    output: wide='tables/segment_match_rates.wide.tsv.gz', long='tables/segment_match_rates.long.tsv.gz', summary='tables/population_match_summary.tsv'
    params: stage='collect_tables', code=CODE
    log: 'logs/collect_tables.log'
    script: TASK
