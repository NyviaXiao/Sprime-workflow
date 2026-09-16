# Draw a two-reference density contour for each configured population/pair.
rule contour:
    input: table='tables/{population}.wide.tsv', script=str(Path(ROOT) / 'workflow/scripts/plot_contour.R')
    output: 'plots/contour/{population}.{pair}.png'
    params: stage='contour', pairs=PAIR_IDS, code=SIG['contour']
    log: 'logs/{population}/contour.{pair}.log'
    script: TASK
