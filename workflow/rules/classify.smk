# Apply the configured Neanderthal and Denisovan thresholds. Every segment not
# assigned to either explicit class is labelled ambiguous.
rule classify:
    input: 'tables/{population}.wide.tsv'
    output: 'classification/{population}.tsv'
    params: stage='classify', settings=C['classification'], code=CODE
    log: 'logs/{population}/classify.log'
    script: TASK
