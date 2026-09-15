# Apply the configured Neanderthal and Denisovan thresholds. Every segment not
# assigned to either explicit class is labelled ambiguous.
rule classify:
    input: 'tables/{population}.wide.tsv'
    output: 'classification/{population}.tsv'
    params: stage='classify', settings=C['classification'], code=SIG['classify']
    log: 'logs/{population}/classify.log'
    script: TASK

# Project-level counts and marker-span coverage for all three classes.
rule classification_summary:
    input:
        classification=expand('classification/{population}.tsv', population=P),
        tables=expand('tables/{population}.wide.tsv', population=P)
    output: 'tables/classification_summary.tsv'
    params: stage='classification_summary', code=SIG['classification_summary']
    log: 'logs/classification_summary.log'
    script: TASK
