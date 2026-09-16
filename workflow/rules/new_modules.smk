"""Independent downstream branches for affinity, adaptive, provenance and report."""

rule affinity_table:
    input: wide='tables/{population}.wide.tsv', classification='classification/{population}.tsv'
    output: 'affinity/{population}/{ref}.tsv.gz'
    params: code=SIG['affinity']
    script: str(Path(ROOT) / 'workflow/scripts/build_affinity_tables.py')

rule landscape:
    input:
        classification='classification/{population}.tsv',
        wide='tables/{population}.wide.tsv',
        aff=lambda w: expand('affinity/{population}/{ref}.tsv.gz', population=w.population, ref=list(R))
    output:
        introgression='plots/landscape/{population}.introgression_landscape.png',
        neanderthal='plots/landscape/{population}.neanderthal_affinity_landscape.png',
        denisovan='plots/landscape/{population}.denisovan_affinity_landscape.png'
    params:
        code=SIG['landscape'], chromosomes=CH, neanderthal_refs=NEAN,
        denisovan_refs=DENI
    script: str(Path(ROOT) / 'workflow/scripts/plot_landscape.R')

rule adaptive_chr:
    input:
        score='sprime/{population}/chr{chrom}.score', vcf='work/{population}/chr{chrom}.vcf.gz',
        samples=C['samples'], mscore=lambda w: expand('archaic_match/{population}/{ref}/chr{chrom}.mscore', population=w.population, ref=list(R), chrom=w.chrom)
    output:
        variants='work/{population}/adaptive.{chrom}.variants.tsv.gz',
        core='work/{population}/adaptive.{chrom}.core.tsv.gz',
        segments='work/{population}/adaptive.{chrom}.segments.tsv'
    params:
        mode='chrom', adaptive=C['adaptive'], references=list(R), tags={k:v['tag'] for k,v in R.items()}, groups={k:v['group'] for k,v in R.items()},
        neanderthal_refs=NEAN, denisovan_refs=DENI, code=SIG['adaptive']
    script: str(Path(ROOT) / 'workflow/scripts/run_adaptive.py')

rule collect_adaptive:
    input: lambda w: expand(['work/{population}/adaptive.{chrom}.variants.tsv.gz','work/{population}/adaptive.{chrom}.core.tsv.gz','work/{population}/adaptive.{chrom}.segments.tsv'], population=P, chrom=CH)
    output:
        variants='adaptive/variants.tsv.gz', core='adaptive/core_variants.tsv.gz', segments='adaptive/segment_summary.tsv', top2='adaptive/top2_candidates.tsv', shared='adaptive/shared_top2_regions.tsv'
    params: mode='collect', outdir='.', top_regions=C['adaptive']['top_regions'], code=SIG['adaptive']
    script: str(Path(ROOT) / 'workflow/scripts/run_adaptive.py')

rule provenance_base:
    input: config='resolved_config.json'
    output:
        run='provenance/run_manifest.json', resolved='provenance/resolved_config.json', software='provenance/software_versions.tsv', inputs='provenance/input_manifest.tsv'
    params: mode='base', root=ROOT, code=SIG['provenance']
    script: str(Path(ROOT) / 'workflow/scripts/build_provenance.py')

rule report:
    input:
        run='provenance/run_manifest.json', classification='tables/classification_summary.tsv', individual='individual_calls/individual_summary.tsv' if C['individual']['enabled'] else [],
        intro=expand('plots/landscape/{population}.introgression_landscape.png', population=P) if C['landscape']['enabled'] else [],
        affinity=expand('plots/landscape/{population}.{kind}.png', population=P, kind=['neanderthal_affinity_landscape','denisovan_affinity_landscape']) if C['landscape']['enabled'] else [],
        contour=lambda w: f"plots/contour/{{P}}.{C['report']['contour_pair'][0]}__{C['report']['contour_pair'][1]}.png".format(P=P[0]),
        gmm='gmm/model_selection.tsv' if C['gmm']['enabled'] else [], gmm_plots=expand('gmm/plots/{population}.{ref}.png', population=P, ref=C['gmm']['target_references']) if C['gmm']['enabled'] else [], adaptive='adaptive/top2_candidates.tsv' if C['adaptive']['enabled'] else []
    output: html='report/report.html', classification_plot='report/classification_summary.png', individual_plot='report/individual_summary.png'
    params: adaptive=C['adaptive']['enabled'], code=SIG['report']
    script: str(Path(ROOT) / 'workflow/scripts/build_report.py')

rule provenance_outputs:
    input: artifacts=FINAL_ARTIFACTS + ['report/report.html']
    output: 'provenance/output_manifest.tsv'
    params: mode='outputs', code=SIG['provenance']
    script: str(Path(ROOT) / 'workflow/scripts/build_provenance.py')
