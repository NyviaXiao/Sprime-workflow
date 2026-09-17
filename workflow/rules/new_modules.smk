"""Independent downstream analysis and presentation branches."""

rule affinity_table:
    input: wide='tables/{population}.wide.tsv', classification='classification/{population}.tsv'
    output: 'affinity/{population}/{ref}.tsv.gz'
    params: code=SIG['affinity_table']
    script: str(Path(ROOT) / 'workflow/scripts/build_affinity_tables.py')

rule introgression_landscape:
    input: classification='classification/{population}.tsv'
    output: introgression='plots/landscape/{population}.introgression_landscape.png'
    params: chromosomes=CH, genome_build=C['genome_build'], code=SIG['introgression_landscape']
    script: str(Path(ROOT) / 'workflow/scripts/plot_landscape.R')

rule affinity_landscape:
    input:
        nean_aff=lambda w: expand('affinity/{population}/{ref}.tsv.gz', population=w.population, ref=NEAN),
        deni_aff=lambda w: expand('affinity/{population}/{ref}.tsv.gz', population=w.population, ref=DENI)
    output:
        neanderthal='plots/landscape/{population}.neanderthal_affinity_landscape.png',
        denisovan='plots/landscape/{population}.denisovan_affinity_landscape.png'
    params:
        chromosomes=CH, genome_build=C['genome_build'], neanderthal_refs=NEAN, denisovan_refs=DENI,
        code=SIG['affinity_landscape']
    script: str(Path(ROOT) / 'workflow/scripts/plot_affinity_landscape.R')

rule adaptive_chr:
    input:
        score='sprime/{population}/chr{chrom}.score', vcf='work/{population}/chr{chrom}.vcf.gz',
        samples=C['samples'],
        mscore=lambda w: expand('archaic_match/{population}/{ref}/chr{chrom}.mscore', population=w.population, ref=list(R), chrom=w.chrom)
    output:
        variants='work/{population}/adaptive.{chrom}.variants.tsv.gz',
        core='work/{population}/adaptive.{chrom}.core.tsv.gz',
        segments='work/{population}/adaptive.{chrom}.segments.tsv'
    params:
        adaptive=C['adaptive'], references=list(R), tags={key:value['tag'] for key,value in R.items()},
        groups={key:value['group'] for key,value in R.items()}, neanderthal_refs=NEAN,
        denisovan_refs=DENI, code=SIG['adaptive_chr']
    script: str(Path(ROOT) / 'workflow/scripts/run_adaptive.py')

rule collect_adaptive:
    input: lambda w: expand(['work/{population}/adaptive.{chrom}.variants.tsv.gz', 'work/{population}/adaptive.{chrom}.core.tsv.gz', 'work/{population}/adaptive.{chrom}.segments.tsv'], population=P, chrom=CH)
    output:
        variants='adaptive/variants.tsv.gz', core='adaptive/core_variants.tsv.gz',
        segments='adaptive/segment_summary.tsv', top2='adaptive/top2_candidates.tsv',
        shared='adaptive/shared_top2_regions.tsv'
    params: code=SIG['adaptive_collect']
    script: str(Path(ROOT) / 'workflow/scripts/collect_adaptive.py')

rule provenance_base:
    input: config='resolved_config.json', validation='qc/validation.tsv', run_info='run_info.json'
    output:
        run='provenance/run_manifest.json', resolved='provenance/resolved_config.json',
        software='provenance/software_versions.tsv', inputs='provenance/input_manifest.tsv'
    params: mode='base', git_commit=CURRENT_GIT_COMMIT, code=SIG['provenance']
    script: str(Path(ROOT) / 'workflow/scripts/build_provenance.py')

rule report:
    input:
        run='provenance/run_manifest.json', software='provenance/software_versions.tsv',
        classification='tables/classification_summary.tsv',
        intro=expand('plots/landscape/{population}.introgression_landscape.png', population=P) if C['landscape']['enabled'] else [],
        affinity=expand('plots/landscape/{population}.{kind}.png', population=P, kind=['neanderthal_affinity_landscape','denisovan_affinity_landscape']) if C['affinity']['enabled'] else [],
        contour=f"plots/contour/{C['report']['contour_population']}.{C['report']['contour_pair'][0]}__{C['report']['contour_pair'][1]}.png" if C['report']['enabled'] and C['plots']['contour']['enabled'] else [],
        individual='individual_calls/individual_summary.tsv' if C['individual']['enabled'] else [],
        gmm='gmm/model_selection.tsv' if C['gmm']['enabled'] else [],
        gmm_plots=expand('gmm/plots/{population}.{ref}.png', population=P, ref=C['gmm']['target_references']) if C['gmm']['enabled'] else [],
        adaptive='adaptive/top2_candidates.tsv' if C['adaptive']['enabled'] else [],
        adaptive_shared='adaptive/shared_top2_regions.tsv' if C['adaptive']['enabled'] else []
    output: ['report/report.html', 'report/report.pdf', 'report/classification_summary.png'] + (['report/individual_summary.png'] if C['individual']['enabled'] else [])
    params: code=SIG['report']
    script: str(Path(ROOT) / 'workflow/scripts/build_report.py')

rule provenance_outputs:
    input: artifacts=FINAL_ARTIFACTS
    output: 'provenance/output_manifest.tsv'
    params: mode='outputs', code=SIG['provenance']
    script: str(Path(ROOT) / 'workflow/scripts/build_provenance.py')
