#!/usr/bin/env python3
"""Render a deterministic HTML/PDF research report from final workflow files.

This presentation layer reads final result tables only. It does not recompute
scientific classifications, GMM decisions, individual calls, or adaptive filters.
"""
import csv
import html
import json
import os
from pathlib import Path
from statistics import mean, median

CLASSES = ('neanderthal', 'denisovan', 'ambiguous')
GMM_COLUMNS = ['population', 'reference', 'n_segments', 'status', 'method',
               'selected_components', 'adjusted_p_1_vs_2', 'second_comparison',
               'adjusted_p_second', 'converged']


def _read_tsv(path):
    if not path:
        return []
    opener = __import__('gzip').open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def _number(value, digits=2):
    try:
        return f'{float(value):.{digits}f}'
    except (TypeError, ValueError):
        return str(value)


def _relative(path, output_dir):
    return html.escape(os.path.relpath(Path(path), output_dir).replace('\\', '/'))


def _table(rows, fields=None):
    if not rows:
        return '<p class="muted">No rows.</p>'
    fields = fields or list(rows[0])
    header = ''.join(f'<th>{html.escape(field)}</th>' for field in fields)
    text = {'population', 'archaic_class', 'adaptive_group', 'chromosome',
            'segment_id', 'reference', 'software', 'version', 'module', 'setting', 'value', 'status', 'method',
            'second_comparison', 'enabled_modules', 'git_commit', 'git_dirty'}
    body = []
    for row in rows:
        cells = [f'<td class="{"text" if field in text else "numeric"}">{html.escape(str(row.get(field, "")))}</td>' for field in fields]
        body.append('<tr>' + ''.join(cells) + '</tr>')
    return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def _figure(title, path, caption, output_dir):
    if not path:
        return ''
    return (f'<figure><h3>{html.escape(title)}</h3><img src="{_relative(path, output_dir)}" '
            f'alt="{html.escape(title)}"><figcaption>{html.escape(caption)}</figcaption></figure>')


def _barplot(rows, value, path, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    populations = sorted({row['population'] for row in rows})
    figure, axis = plt.subplots(figsize=(8.4, 4.4))
    palette = {'neanderthal': '#3f6f9f', 'denisovan': '#a95243', 'ambiguous': '#7d8790'}
    for index, archaic_class in enumerate(CLASSES):
        values = [float(next((row[value] for row in rows if row['population'] == pop and row['archaic_class'] == archaic_class), 0)) for pop in populations]
        axis.bar([item + (index - 1) * .23 for item in range(len(populations))], values,
                 .23, label=archaic_class.capitalize(), color=palette[archaic_class])
    axis.set_xticks(range(len(populations)), populations)
    axis.set_ylabel(value.replace('_', ' ')); axis.set_title(title, weight='semibold')
    axis.spines[['top', 'right']].set_visible(False); axis.legend(frameon=False, ncol=3)
    figure.tight_layout(); figure.savefig(path, dpi=200); plt.close(figure)


def _classification_narrative(rows):
    result = []
    for population in sorted({row['population'] for row in rows}):
        values = {klass: next((row for row in rows if row['population'] == population and row['archaic_class'] == klass), {}) for klass in CLASSES}
        counts = [values[klass].get('n_segments', 0) for klass in CLASSES]
        coverage = [_number(values[klass].get('genome_coverage_mb', 0)) for klass in CLASSES]
        result.append(f'{population} contained {counts[0]} Neanderthal-classified, {counts[1]} Denisovan-classified, and {counts[2]} ambiguous SPrime segments; corresponding marker-supported genomic coverage was {coverage[0]}, {coverage[1]}, and {coverage[2]} Mb.')
    return ' '.join(result)


def _individual_narrative(rows):
    if not rows:
        return 'Not enabled for this run.'
    sentences = []
    for row in rows:
        sentences.append(f"{row['population']} showed {_number(row['per_individual_introgressed_mb'])} Mb of marker-supported {row['archaic_class']}-classified haplotype sequence per individual on average, with {row['n_introgressed_individuals']} individuals carrying at least one such tract.")
    return ' '.join(sentences) + ' These are marker-supported haplotype tracts and are not inferred recombination breakpoints.'


def _gmm_narrative(rows):
    if not rows:
        return 'Not enabled for this run.'
    selected = [row.get('selected_components', '') for row in rows]
    counts = {item: sum(value == str(item) for value in selected) for item in (1, 2, 3)}
    insufficient = sum(value in ('', 'NA', 'None') for value in selected)
    return (f'Of {len(rows)} population-reference models, {counts[1]} selected one component, {counts[2]} selected two components, and {counts[3]} selected three components; {insufficient} analyses had insufficient data. Mixture component counts describe statistical structure in the match-rate distribution and are not interpreted here as direct counts of historical introgression events.')


def _adaptive_narrative(segment_rows, top2_rows):
    if not segment_rows:
        return 'Not enabled for this run.'
    nean = sum(int(row.get('neanderthal_pass', 0)) == 1 for row in segment_rows)
    deni = sum(int(row.get('denisovan_pass', 0)) == 1 for row in segment_rows)
    sentence = (f'Adaptive screening retained {nean} Neanderthal-associated and {deni} Denisovan-associated passing segments. Candidates were ranked independently within each archaic group using core mean introgressed-allele frequency. The report displays up to two candidates per group and population.')
    short = []
    for population in sorted({row['population'] for row in top2_rows}):
        for group in ('neanderthal', 'denisovan'):
            if len([row for row in top2_rows if row['population'] == population and row.get('adaptive_group') == group]) == 1:
                short.append(f'Only one qualifying candidate was available for {population} {group}.')
    return sentence + (' ' + ' '.join(short) if short else '')


def _affinity_narrative(paths, cfg):
    if not paths:
        return 'Not enabled for this run.'
    refs = {ref['id']: ref for ref in cfg['archaic_references']}
    summaries = []
    for path in paths:
        ref_id = Path(path).stem.replace('.tsv', '')
        ref = refs.get(ref_id)
        if not ref:
            continue
        rows = [row for row in _read_tsv(path) if row.get('archaic_class') == ref['group']]
        rates = []
        for row in rows:
            try:
                rates.append(float(row['match_rate']))
            except (KeyError, TypeError, ValueError):
                pass
        population = rows[0]['population'] if rows else Path(path).parent.name
        if rates:
            summaries.append((population, ref['group'], ref['tag'], len(rows), median(rates), mean(rates)))
    if not summaries:
        return 'No finite reference-specific match rates were available among the corresponding classified segments.'
    sentences = []
    for population in sorted({item[0] for item in summaries}):
        for group in ('neanderthal', 'denisovan'):
            items = [item for item in summaries if item[0] == population and item[1] == group]
            if items:
                details = ', '.join(f'{tag} (n={count}; median={_number(med)}; mean={_number(avg)})' for _, _, tag, count, med, avg in items)
                sentences.append(f'Among {population} {group}-classified segments, reference-specific match-rate summaries were {details}.')
    return ' '.join(sentences)


def _settings_rows(cfg):
    rows = [('SPrime', 'minimum score', cfg['sprime_minscore']), ('Summary', 'minimum callable', cfg['summary_min_callable'])]
    c = cfg['classification']
    rows.extend([('Classification', 'Neanderthal references', ', '.join(c['neanderthal_references'])),
                 ('Classification', 'Denisovan references', ', '.join(c['denisovan_references'])),
                 ('Classification', 'Neanderthal threshold / Denisovan exclusion', f"> {c['neanderthal_gt']} / < {c['denisovan_lt_for_neanderthal']}"),
                 ('Classification', 'Denisovan threshold / Neanderthal exclusion', f"> {c['denisovan_gt']} / < {c['neanderthal_lt_for_denisovan']}")])
    if cfg['gmm']['enabled']:
        g = cfg['gmm']; rows.extend([('GMM', 'score_gt / min_callable / min_segments', f"{g['score_gt']} / {g['min_callable']} / {g['min_segments']}"), ('GMM', 'method / alpha', f"{g['method']} / {g['alpha']}"), ('GMM', 'target references', ', '.join(g['target_references']))])
    if cfg['individual']['enabled']:
        rows.append(('Individual', 'minimum markers', cfg['individual']['min_markers']))
    if cfg['adaptive']['enabled']:
        a = cfg['adaptive']; rows.extend([('Adaptive', 'minimum allele frequency / maximum frequency drop', f"{a['min_allele_frequency']} / {a['max_frequency_drop']}"), ('Adaptive', 'minimum callable sites', a['min_callable_sites']), ('Adaptive', 'Neanderthal / Denisovan match thresholds', f"{a['neanderthal_match_threshold']} / {a['denisovan_match_threshold']}")])
    return [{'module': module, 'setting': setting, 'value': value} for module, setting, value in rows]


def _selected_gmm_plots(rows, plots):
    chosen = {f"{row.get('population')}__{row.get('reference')}" for row in rows if row.get('selected_components', '') in ('2', '3')}
    return [path for path in plots if Path(path).stem.replace('.', '__') in chosen]


def _write_pdf(html_path, pdf_path):
    try:
        from weasyprint import HTML
    except ImportError as error:
        raise RuntimeError('Report PDF generation requires WeasyPrint; install the environment.yaml dependency.') from error
    HTML(filename=str(html_path), base_url=str(html_path.parent.resolve())).write_pdf(str(pdf_path))


def build_report(output_html, output_pdf, run_manifest, software_versions, resolved_config,
                 classification, intro_images, affinity_images, affinity_tables, contour,
                 individual, gmm, gmm_plots, adaptive_segments, adaptive_top2,
                 adaptive_shared=None, classification_plot=None, individual_plot=None):
    """Write formal HTML/PDF report output from final result artifacts."""
    output_html, output_pdf = Path(output_html), Path(output_pdf)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(Path(run_manifest).read_text(encoding='utf-8'))
    cfg = json.loads(Path(resolved_config).read_text(encoding='utf-8'))
    software, class_rows = _read_tsv(software_versions), _read_tsv(classification)
    individual_rows, gmm_rows = _read_tsv(individual), _read_tsv(gmm)
    segment_rows, top2_rows, shared_rows = _read_tsv(adaptive_segments), _read_tsv(adaptive_top2), _read_tsv(adaptive_shared)
    if classification_plot:
        _barplot(class_rows, 'genome_coverage_mb', classification_plot, 'Classification genome coverage (Mb)')
    if individual_rows and individual_plot:
        _barplot(individual_rows, 'per_individual_introgressed_mb', individual_plot, 'Per-individual archaic haplotype amount (Mb)')
    outdir, enabled = output_html.parent, manifest.get('enabled_modules', {})
    cards = ''.join(f'<div class="card"><span>{html.escape(label)}</span><strong>{html.escape(str(value))}</strong></div>' for label, value in (
        ('Genome build', manifest.get('genome_build', '')), ('Populations', ', '.join(manifest.get('populations', []))), ('Archaic references', len(manifest.get('archaic_references', []))), ('Workflow commit', manifest.get('git_commit', '')), ('Working tree dirty', manifest.get('git_dirty', False)), ('Generated', manifest.get('run_time', ''))))
    executive = '<p>' + _classification_narrative(class_rows) + '</p>'
    if individual_rows: executive += '<p>' + _individual_narrative(individual_rows) + '</p>'
    if gmm_rows: executive += f'<p>GMM analysis was performed for {len(gmm_rows)} population-reference models; {sum(row.get("selected_components", "") in ("2", "3") for row in gmm_rows)} selected more than one mixture component.</p>'
    if segment_rows: executive += '<p>' + _adaptive_narrative(segment_rows, top2_rows) + '</p>'
    intro_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} - Introgression Landscape', path, 'Figure C1. Genome-wide distribution of classified SPrime segments. Colored intervals are overlaid directly on GRCh37 chromosome bodies; grey regions are not covered by displayed classified segments.', outdir) for path in intro_images) or '<p class="muted">Not enabled for this run.</p>'
    nean_paths = [path for path in affinity_images if '.neanderthal_affinity_landscape.' in str(path)]
    deni_paths = [path for path in affinity_images if '.denisovan_affinity_landscape.' in str(path)]
    nean_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} - Neanderthal Affinity Landscape', path, 'Figure D1. Reference-specific affinity of Neanderthal-classified segments. Each chromosome body is internally divided into reference layers; color intensity represents segment match rate from 0 to 1.', outdir) for path in nean_paths) or '<p class="muted">Not enabled for this run.</p>'
    deni_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} - Denisovan Affinity Landscape', path, 'Figure D2. Reference-specific affinity of Denisovan-classified segments. Each chromosome body is internally divided into reference layers; color intensity represents segment match rate from 0 to 1.', outdir) for path in deni_paths) or '<p class="muted">Not enabled for this run.</p>'
    compact_adaptive = ['population', 'adaptive_group', 'adaptive_rank', 'chromosome', 'segment_id', 'candidate_start', 'candidate_end', 'core_mean_AF', 'max_AF', 'core_variant_count', 'neanderthal_pass', 'denisovan_pass']
    compact_gmm = [field for field in GMM_COLUMNS if gmm_rows and field in gmm_rows[0]]
    refs = {ref['id']: ref['tag'] for ref in cfg['archaic_references']}
    contour_pair = cfg['report'].get('contour_pair', [])
    contour_text = (f"Figure E1 compares segment-level match rates to {refs.get(contour_pair[0], contour_pair[0])} and {refs.get(contour_pair[1], contour_pair[1])} in {cfg['report'].get('contour_population', '')}." if contour else 'Not enabled for this run.')
    sections = [
        ('overview', 'A. Run Overview', '<p>Run metadata and software versions identify this analysis without embedding filesystem paths or the complete resolved configuration.</p>' + _table([{'run_time': manifest.get('run_time', ''), 'genome_build': manifest.get('genome_build', ''), 'populations': ', '.join(manifest.get('populations', [])), 'chromosomes': ', '.join(manifest.get('chromosomes', [])), 'archaic_references': ', '.join(manifest.get('archaic_references', [])), 'enabled_modules': json.dumps(enabled, sort_keys=True), 'git_commit': manifest.get('git_commit', ''), 'git_dirty': manifest.get('git_dirty', False)}]) + '<h3>A2. Analysis Settings</h3>' + _table(_settings_rows(cfg)) + '<h3>Software environment</h3>' + _table(software)),
        ('classification', 'B. Classification Results', '<p>' + _classification_narrative(class_rows) + '</p>' + _table(class_rows) + (_figure('Classification genome coverage', classification_plot, 'Figure B1. Marker-supported genomic coverage by final segment classification.', outdir) if classification_plot else '')),
        ('landscape', 'C. Introgression Landscape', '<p>Classified segments are shown at their genomic positions on a common GRCh37 chromosome background; this display does not infer hotspots or enrichment.</p>' + intro_html),
        ('affinity', 'D. Archaic Affinity Analysis', '<p>' + _affinity_narrative(affinity_tables, cfg) + '</p><p>Affinity plots are descriptive reference-specific match rates for segments in the corresponding final class. They do not assign a donor or source population.</p><h3>D1. Neanderthal Affinity</h3>' + nean_html + '<h3>D2. Denisovan Affinity</h3>' + deni_html),
        ('contour', 'E. Pairwise Archaic Affinity Contour', '<p>' + contour_text + '</p>' + (_figure('Configured Neanderthal x Denisovan contour', contour, 'Figure E1. Segment-level pairwise match-rate density. The plot is descriptive and does not assign archaic source identity.', outdir) if contour else '')),
        ('individual', 'F. Individual-level Introgression', ('<p>' + _individual_narrative(individual_rows) + '</p>' + _table(individual_rows) + (_figure('Per-individual archaic haplotype amount', individual_plot, 'Figure F1. Mean marker-supported archaic haplotype sequence amount per target individual.', outdir) if individual_plot else '')) if individual_rows else '<p class="muted">Not enabled for this run.</p>'),
        ('gmm', 'G. GMM Analysis', '<p>' + _gmm_narrative(gmm_rows) + '</p>' + (_table(gmm_rows, compact_gmm) + ''.join(_figure(Path(path).stem, path, 'Figure G. GMM result shown only where more than one component was selected.', outdir) for path in _selected_gmm_plots(gmm_rows, gmm_plots)) if gmm_rows else '')),
        ('adaptive', 'H. Adaptive Introgression Candidates', '<p>' + _adaptive_narrative(segment_rows, top2_rows) + '</p><h3>H1. Neanderthal-associated candidates</h3>' + (_table([row for row in top2_rows if row.get('adaptive_group') == 'neanderthal'], compact_adaptive) if any(row.get('adaptive_group') == 'neanderthal' for row in top2_rows) else '<p class="muted">No qualifying candidate.</p>') + '<h3>H2. Denisovan-associated candidates</h3>' + (_table([row for row in top2_rows if row.get('adaptive_group') == 'denisovan'], compact_adaptive) if any(row.get('adaptive_group') == 'denisovan' for row in top2_rows) else '<p class="muted">No qualifying candidate.</p>') + '<h3>H3. Shared Top2 regions</h3>' + (_table(shared_rows) if shared_rows else '<p class="muted">No shared Top2 regions were identified.</p>')),
        ('provenance', 'I. Reproducibility & Provenance', '<p>Detailed resolved configuration and input/output manifests remain in the provenance directory.</p>' + _table([{'git_commit': manifest.get('git_commit', ''), 'git_dirty': manifest.get('git_dirty', False), 'git_diff_sha256': manifest.get('git_diff_sha256', ''), 'genome_build': manifest.get('genome_build', ''), 'populations': ', '.join(manifest.get('populations', [])), 'chromosomes': ', '.join(manifest.get('chromosomes', [])), 'archaic_references': ', '.join(manifest.get('archaic_references', [])), 'enabled_modules': json.dumps(enabled, sort_keys=True)}]) + '<h3>Core tool versions</h3>' + _table([row for row in software if row.get('software') in ('snakemake', 'bcftools', 'java', 'R')])),
    ]
    toc = ''.join(f'<a href="#{anchor}">{title}</a>' for anchor, title, _ in sections)
    rendered = ''.join(f'<section id="{anchor}"><h2>{title}</h2>{content}</section>' for anchor, title, content in sections)
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>SPrime Archaic Introgression Analysis Report</title><style>
@page {{ size: A4; margin: 16mm 15mm 17mm 15mm; }} * {{ box-sizing:border-box; }} body {{ margin:0; background:#fff; color:#24313a; font-family:Arial,Helvetica,sans-serif; font-size:14px; line-height:1.55; }} .document {{ max-width:1320px; margin:0 auto; padding:42px 52px 64px; }} .hero {{ border-bottom:2px solid #dbe3e8; padding-bottom:24px; }} h1 {{ font-size:30px; margin:0 0 4px; color:#18252d; }} .subtitle {{ color:#60717d; margin:0 0 22px; }} .cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; }} .card {{ background:#f4f7f8; border:1px solid #e0e7eb; border-radius:7px; padding:12px; min-height:72px; break-inside:avoid; page-break-inside:avoid; }} .card span {{ display:block; color:#637580; font-size:11px; text-transform:uppercase; letter-spacing:.04em; }} .card strong {{ display:block; margin-top:4px; overflow-wrap:anywhere; }} .toc {{ margin:26px 0 36px; padding:16px 20px; background:#f7f9fa; border-left:4px solid #718c9e; }} .toc h2 {{ margin-top:0; font-size:16px; }} .toc a {{ display:inline-block; margin:3px 16px 3px 0; color:#365b71; text-decoration:none; }} section {{ border-top:1px solid #dbe3e8; margin-top:34px; padding-top:21px; }} h2 {{ font-size:21px; color:#1c303c; margin:0 0 12px; }} h3 {{ font-size:16px; color:#365467; margin:22px 0 9px; }} .table-wrap {{ overflow-x:auto; margin:12px 0 20px; }} table {{ width:100%; border-collapse:collapse; background:white; }} th {{ background:#eaf0f3; color:#29414f; font-weight:600; }} th,td {{ padding:7px 8px; border:1px solid #d8e0e4; vertical-align:top; }} tbody tr:nth-child(even) {{ background:#fafcfc; }} td.numeric {{ text-align:right; font-variant-numeric:tabular-nums; }} td.text {{ text-align:left; }} figure {{ margin:22px auto 30px; padding:12px 0; break-inside:avoid; page-break-inside:avoid; }} figure h3 {{ margin-top:0; }} figure img {{ display:block; max-width:100%; max-height:245mm; margin:0 auto; }} figcaption {{ color:#5c6a73; font-size:12px; margin:9px auto 0; max-width:980px; }} .muted {{ color:#687984; font-style:italic; }} @media print {{ .document {{ max-width:none; padding:0; }} .hero,.toc,figure,.card {{ break-inside:avoid; page-break-inside:avoid; }} h2 {{ break-after:avoid; }} }}</style></head><body><main class="document"><header class="hero"><h1>SPrime Archaic Introgression Analysis Report</h1><p class="subtitle">Final workflow results and reproducibility record</p><div class="cards">{cards}</div></header><section id="executive"><h2>Executive Summary</h2>{executive}</section><nav class="toc"><h2>Table of Contents</h2>{toc}</nav>{rendered}</main></body></html>'''
    output_html.write_text(page, encoding='utf-8'); _write_pdf(output_html, output_pdf)


if 'snakemake' in globals():
    outputs = list(snakemake.output)
    build_report(outputs[0], outputs[1], snakemake.input.run, snakemake.input.software,
                 snakemake.input.resolved_config, snakemake.input.classification,
                 list(snakemake.input.intro), list(snakemake.input.affinity),
                 list(snakemake.input.affinity_tables), snakemake.input.contour or None,
                 snakemake.input.individual or None, snakemake.input.gmm or None,
                 list(snakemake.input.gmm_plots), snakemake.input.adaptive_segments or None,
                 snakemake.input.adaptive_top2 or None, snakemake.input.adaptive_shared or None,
                 outputs[2], outputs[3] if len(outputs) > 3 else None)
