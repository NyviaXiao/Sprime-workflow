#!/usr/bin/env python3
"""Render a deterministic HTML/PDF research report from final workflow tables.

This module never recalculates biological quantities. Its prose is a compact
rendering of final tables and avoids donor, event-count, or mechanism claims.
"""
import csv
import html
import json
import os
from pathlib import Path


CLASSES = ('neanderthal', 'denisovan', 'ambiguous')


def _read_tsv(path):
    if not path:
        return []
    with open(path, newline='', encoding='utf-8') as handle:
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
    body = []
    text_fields = {'population', 'archaic_class', 'adaptive_group', 'chromosome',
                   'segment_id', 'reference', 'software', 'version'}
    for row in rows:
        cells = []
        for field in fields:
            kind = 'text' if field in text_fields else 'numeric'
            cells.append(f'<td class="{kind}">{html.escape(str(row.get(field, "")))}</td>')
        body.append('<tr>' + ''.join(cells) + '</tr>')
    return f'<div class="table-wrap"><table><thead><tr>{header}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def _figure(title, path, caption, output_dir):
    if not path:
        return ''
    return (f'<figure><h3>{html.escape(title)}</h3>'
            f'<img src="{_relative(path, output_dir)}" alt="{html.escape(title)}">'
            f'<figcaption>{html.escape(caption)}</figcaption></figure>')


def _barplot(rows, value, path, title):
    """Create compact report-local summary plots from existing final tables."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    populations = sorted({row['population'] for row in rows})
    figure, axis = plt.subplots(figsize=(8.4, 4.4))
    palette = {'neanderthal': '#3f6f9f', 'denisovan': '#a95243', 'ambiguous': '#7d8790'}
    width = .23
    for index, archaic_class in enumerate(CLASSES):
        values = [float(next((row[value] for row in rows if row['population'] == pop and row['archaic_class'] == archaic_class), 0)) for pop in populations]
        axis.bar([item + (index - 1) * width for item in range(len(populations))], values,
                 width, label=archaic_class.capitalize(), color=palette[archaic_class])
    axis.set_xticks(range(len(populations)), populations)
    axis.set_ylabel(value.replace('_', ' '))
    axis.set_title(title, weight='semibold')
    axis.spines[['top', 'right']].set_visible(False)
    axis.legend(frameon=False, ncol=3)
    figure.tight_layout()
    figure.savefig(path, dpi=200)
    plt.close(figure)


def _classification_narrative(rows):
    sentences = []
    for population in sorted({row['population'] for row in rows}):
        values = {klass: next((row for row in rows if row['population'] == population and row['archaic_class'] == klass), {}) for klass in CLASSES}
        counts = [str(values[klass].get('n_segments', 0)) for klass in CLASSES]
        coverage = [_number(values[klass].get('genome_coverage_mb', 0), 2) for klass in CLASSES]
        sentences.append(f'{population} contained {counts[0]} Neanderthal-classified, {counts[1]} Denisovan-classified, and {counts[2]} ambiguous SPrime segments; corresponding marker-supported genomic coverage was {coverage[0]}, {coverage[1]}, and {coverage[2]} Mb.')
    return ' '.join(sentences)


def _gmm_narrative(rows):
    if not rows:
        return 'Not enabled for this run.'
    selected = [row.get('selected_components', '') for row in rows]
    counts = {component: sum(value == str(component) for value in selected) for component in (1, 2, 3)}
    insufficient = sum(value in ('', 'NA', 'None') for value in selected)
    return (f'Of {len(rows)} population-reference models, {counts[1]} selected one component, {counts[2]} selected two components, and {counts[3]} selected three components; {insufficient} analyses had insufficient data. Mixture component counts describe statistical structure in the match-rate distribution and are not interpreted here as direct counts of historical introgression events.')


def _adaptive_narrative(rows):
    if not rows:
        return 'Not enabled for this run.'
    counts = {group: sum(row.get('adaptive_group') == group for row in rows) for group in ('neanderthal', 'denisovan')}
    return (f'After allele-frequency and archaic-matching filters, {counts["neanderthal"]} Neanderthal-associated and {counts["denisovan"]} Denisovan-associated candidate rows were retained for independent group-specific Top2 ranking. A group with fewer than two qualifying segments is reported without padding.')


def _selected_gmm_plots(rows, plots):
    chosen = {f"{row.get('population')}__{row.get('reference')}" for row in rows
              if row.get('selected_components', '') not in ('', 'NA', 'None') and int(row['selected_components']) > 1}
    return [path for path in plots if Path(path).stem.replace('.', '__') in chosen]


def _write_pdf(html_path, pdf_path):
    try:
        from weasyprint import HTML
    except ImportError as error:
        raise RuntimeError('Report PDF generation requires WeasyPrint; install the environment.yaml dependency.') from error
    HTML(filename=str(html_path), base_url=str(html_path.parent.resolve())).write_pdf(str(pdf_path))


def build_report(output_html, output_pdf, run_manifest, software_versions, classification,
                 intro_images, affinity_images, contour, individual, gmm, gmm_plots,
                 adaptive, adaptive_shared=None, classification_plot=None, individual_plot=None):
    """Write a formal report and PDF without performing scientific computation."""
    output_html, output_pdf = Path(output_html), Path(output_pdf)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(Path(run_manifest).read_text(encoding='utf-8'))
    software, class_rows = _read_tsv(software_versions), _read_tsv(classification)
    if classification_plot:
        _barplot(class_rows, 'genome_coverage_mb', classification_plot, 'Classification genome coverage (Mb)')
    individual_rows = _read_tsv(individual)
    if individual_rows and individual_plot:
        _barplot(individual_rows, 'per_individual_introgressed_mb', individual_plot, 'Per-individual archaic haplotype amount (Mb)')
    gmm_rows, adaptive_rows, shared_rows = _read_tsv(gmm), _read_tsv(adaptive), _read_tsv(adaptive_shared)
    outdir, enabled = output_html.parent, manifest.get('enabled_modules', {})
    cards = ''.join(f'<div class="card"><span>{html.escape(label)}</span><strong>{html.escape(str(value))}</strong></div>' for label, value in (
        ('Genome build', manifest.get('genome_build', '')), ('Populations', ', '.join(manifest.get('populations', []))),
        ('Archaic references', len(manifest.get('archaic_references', []))), ('Workflow commit', manifest.get('git_commit', '')),
        ('Generated', manifest.get('run_time', ''))))
    executive = '<p>' + _classification_narrative(class_rows) + '</p>'
    if individual_rows:
        executive += '<p>Individual-level calling produced marker-supported archaic haplotype tract summaries for the enabled populations.</p>'
    if gmm_rows:
        executive += f'<p>GMM analysis was performed for {len(gmm_rows)} population-reference models; {sum(row.get("selected_components", "") in ("2", "3") for row in gmm_rows)} selected more than one mixture component.</p>'
    if adaptive_rows:
        executive += '<p>' + _adaptive_narrative(adaptive_rows) + '</p>'

    intro_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} — Introgression Landscape', path, 'Figure C1. Genome-wide distribution of classified SPrime segments. Colored intervals are overlaid directly on GRCh37 chromosome bodies; grey regions are not covered by displayed classified segments.', outdir) for path in intro_images) or '<p class="muted">Not enabled for this run.</p>'
    nean_paths = [path for path in affinity_images if '.neanderthal_affinity_landscape.' in str(path)]
    deni_paths = [path for path in affinity_images if '.denisovan_affinity_landscape.' in str(path)]
    nean_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} — Neanderthal Affinity Landscape', path, 'Figure D1. Reference-specific affinity of Neanderthal-classified segments. Each chromosome body is internally divided into reference layers; color intensity represents segment match rate from 0 to 1.', outdir) for path in nean_paths) or '<p class="muted">Not enabled for this run.</p>'
    deni_html = ''.join(_figure(f'{Path(path).stem.split(".")[0]} — Denisovan Affinity Landscape', path, 'Figure D2. Reference-specific affinity of Denisovan-classified segments. Each chromosome body is internally divided into reference layers; color intensity represents segment match rate from 0 to 1.', outdir) for path in deni_paths) or '<p class="muted">Not enabled for this run.</p>'
    compact = ['population', 'adaptive_group', 'adaptive_rank', 'chromosome', 'segment_id', 'candidate_start', 'candidate_end', 'core_mean_AF', 'max_AF', 'core_variant_count', 'neanderthal_pass', 'denisovan_pass']
    nean_candidates = [row for row in adaptive_rows if row.get('adaptive_group') == 'neanderthal']
    deni_candidates = [row for row in adaptive_rows if row.get('adaptive_group') == 'denisovan']
    sections = [
        ('overview', 'A. Run Overview', '<p>Run metadata and software versions identify this analysis without embedding the complete resolved configuration.</p>' + _table([{'run_time': manifest.get('run_time', ''), 'genome_build': manifest.get('genome_build', ''), 'populations': ', '.join(manifest.get('populations', [])), 'chromosomes': ', '.join(manifest.get('chromosomes', [])), 'archaic_references': ', '.join(manifest.get('archaic_references', [])), 'enabled_modules': json.dumps(enabled, sort_keys=True), 'git_commit': manifest.get('git_commit', '')}]) + '<h3>Software environment</h3>' + _table(software)),
        ('classification', 'B. Classification Results', '<p>' + _classification_narrative(class_rows) + '</p>' + _table(class_rows) + (_figure('Classification genome coverage', classification_plot, 'Figure B1. Marker-supported genomic coverage by final segment classification.', outdir) if classification_plot else '')),
        ('landscape', 'C. Introgression Landscape', '<p>Classified segments are shown at their genomic positions on a common GRCh37 chromosome background; this display does not infer hotspots or enrichment.</p>' + intro_html),
        ('affinity', 'D. Archaic Affinity Analysis', '<p>Affinity plots display reference-specific match rates only for segments in the corresponding final class. They do not assign a donor or source population.</p><h3>D1. Neanderthal Affinity</h3>' + nean_html + '<h3>D2. Denisovan Affinity</h3>' + deni_html),
        ('contour', 'E. Pairwise Archaic Affinity Contour', '<p>The contour plot compares segment-level match rates between the configured Neanderthal and Denisovan references for the selected report population.</p>' + (_figure('Configured Neanderthal × Denisovan contour', contour, 'Figure E1. Segment-level pairwise match-rate density. The plot is descriptive and does not assign archaic source identity.', outdir) if contour else '<p class="muted">Not enabled for this run.</p>')),
        ('individual', 'F. Individual-level Introgression', ('<p>The mean marker-supported archaic haplotype amount and the number of target individuals with at least one tract are summarized below. Tracts are haplotype-specific and are not recombination breakpoints.</p>' + _table(individual_rows) + (_figure('Per-individual archaic haplotype amount', individual_plot, 'Figure F1. Mean marker-supported archaic haplotype sequence amount per target individual.', outdir) if individual_plot else '')) if individual_rows else '<p class="muted">Not enabled for this run.</p>'),
        ('gmm', 'G. GMM Analysis', '<p>' + _gmm_narrative(gmm_rows) + '</p>' + (_table(gmm_rows) + ''.join(_figure(Path(path).stem, path, 'Figure G. GMM result shown only where more than one component was selected.', outdir) for path in _selected_gmm_plots(gmm_rows, gmm_plots)) if gmm_rows else '')),
        ('adaptive', 'H. Adaptive Introgression Candidates', '<p>' + _adaptive_narrative(adaptive_rows) + '</p><h3>H1. Neanderthal-associated candidates</h3>' + (_table(nean_candidates, compact) if nean_candidates else '<p class="muted">No qualifying candidate.</p>') + '<h3>H2. Denisovan-associated candidates</h3>' + (_table(deni_candidates, compact) if deni_candidates else '<p class="muted">No qualifying candidate.</p>') + '<h3>H3. Shared Top2 regions</h3>' + (_table(shared_rows) if shared_rows else '<p class="muted">No shared Top2 regions were identified.</p>')),
        ('provenance', 'I. Reproducibility & Provenance', '<p>The full resolved configuration and input/output manifests are stored in the provenance directory. The summary below records the run identity and software environment used for this report.</p>' + _table([{'git_commit': manifest.get('git_commit', ''), 'genome_build': manifest.get('genome_build', ''), 'populations': ', '.join(manifest.get('populations', [])), 'chromosomes': ', '.join(manifest.get('chromosomes', [])), 'archaic_references': ', '.join(manifest.get('archaic_references', [])), 'enabled_modules': json.dumps(enabled, sort_keys=True)}]) + _table(software)),
    ]
    toc = ''.join(f'<a href="#{anchor}">{title}</a>' for anchor, title, _ in sections)
    rendered = ''.join(f'<section id="{anchor}"><h2>{title}</h2>{content}</section>' for anchor, title, content in sections)
    page = f'''<!doctype html><html><head><meta charset="utf-8"><title>SPrime Archaic Introgression Analysis Report</title><style>
@page {{ size: A4; margin: 16mm 15mm 17mm 15mm; }} * {{ box-sizing:border-box; }} body {{ margin:0; background:#fff; color:#24313a; font-family:Arial,Helvetica,sans-serif; font-size:14px; line-height:1.55; }} .document {{ max-width:1320px; margin:0 auto; padding:42px 52px 64px; }} .hero {{ border-bottom:2px solid #dbe3e8; padding-bottom:24px; }} h1 {{ font-size:30px; margin:0 0 4px; color:#18252d; }} .subtitle {{ color:#60717d; margin:0 0 22px; }} .cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; }} .card {{ background:#f4f7f8; border:1px solid #e0e7eb; border-radius:7px; padding:12px; min-height:72px; }} .card span {{ display:block; color:#637580; font-size:11px; text-transform:uppercase; letter-spacing:.04em; }} .card strong {{ display:block; margin-top:4px; overflow-wrap:anywhere; }} .toc {{ margin:26px 0 36px; padding:16px 20px; background:#f7f9fa; border-left:4px solid #718c9e; }} .toc h2 {{ margin-top:0; font-size:16px; }} .toc a {{ display:inline-block; margin:3px 16px 3px 0; color:#365b71; text-decoration:none; }} section {{ border-top:1px solid #dbe3e8; margin-top:34px; padding-top:21px; }} h2 {{ font-size:21px; color:#1c303c; margin:0 0 12px; }} h3 {{ font-size:16px; color:#365467; margin:22px 0 9px; }} .table-wrap {{ overflow-x:auto; margin:12px 0 20px; }} table {{ width:100%; border-collapse:collapse; background:white; }} th {{ background:#eaf0f3; color:#29414f; font-weight:600; }} th,td {{ padding:7px 8px; border:1px solid #d8e0e4; vertical-align:top; }} tbody tr:nth-child(even) {{ background:#fafcfc; }} td.numeric {{ text-align:right; font-variant-numeric:tabular-nums; }} td.text {{ text-align:left; }} figure {{ margin:22px auto 30px; padding:12px 0; break-inside:avoid; page-break-inside:avoid; }} figure h3 {{ margin-top:0; }} figure img {{ display:block; max-width:100%; max-height:245mm; margin:0 auto; }} figcaption {{ color:#5c6a73; font-size:12px; margin:9px auto 0; max-width:980px; }} .muted {{ color:#687984; font-style:italic; }} @media print {{ .document {{ max-width:none; padding:0; }} .hero,.toc,section,figure,.table-wrap {{ break-inside:avoid; page-break-inside:avoid; }} h2 {{ break-after:avoid; }} }}</style></head><body><main class="document"><header class="hero"><h1>SPrime Archaic Introgression Analysis Report</h1><p class="subtitle">Final workflow results and reproducibility record</p><div class="cards">{cards}</div></header><section id="executive"><h2>Executive Summary</h2>{executive}</section><nav class="toc"><h2>Table of Contents</h2>{toc}</nav>{rendered}</main></body></html>'''
    output_html.write_text(page, encoding='utf-8')
    _write_pdf(output_html, output_pdf)


if 'snakemake' in globals():
    outputs = list(snakemake.output)
    build_report(outputs[0], outputs[1], snakemake.input.run, snakemake.input.software, snakemake.input.classification, list(snakemake.input.intro), list(snakemake.input.affinity), snakemake.input.contour or None, snakemake.input.individual or None, snakemake.input.gmm or None, list(snakemake.input.gmm_plots), snakemake.input.adaptive or None, snakemake.input.adaptive_shared or None, outputs[2], outputs[3] if len(outputs) > 3 else None)
