#!/usr/bin/env python3
"""Render the final HTML report from already computed final artifacts."""
import csv
import html
import json
from pathlib import Path

SECTION_ORDER = ['Run overview', 'Classification summary', 'Introgression landscape',
                 'Landscape + affinity', 'Neanderthal × Denisovan contour plot',
                 'Individual summary', 'GMM summary', 'Adaptive candidate summary']


def _read_tsv(path):
    if not path:
        return []
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def _table(rows):
    if not rows:
        return '<p>No rows.</p>'
    fields = list(rows[0])
    header = ''.join(f'<th>{html.escape(field)}</th>' for field in fields)
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(row.get(field, "")))}</td>' for field in fields) + '</tr>' for row in rows)
    return f'<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>'


def _image(path, width=800):
    return f'<img src="../{html.escape(str(path))}" width="{width}">'


def _barplot(rows, value, path, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    populations = sorted({row['population'] for row in rows})
    classes = ['neanderthal', 'denisovan', 'ambiguous']
    figure, axis = plt.subplots(figsize=(8, 4))
    width = .24
    for i, archaic_class in enumerate(classes):
        values = [float(next((row[value] for row in rows if row['population'] == pop and row['archaic_class'] == archaic_class), 0)) for pop in populations]
        axis.bar([index + (i - 1) * width for index in range(len(populations))], values,
                 width, label=archaic_class)
    axis.set_xticks(range(len(populations))); axis.set_xticklabels(populations)
    axis.set_ylabel(value); axis.set_title(title); axis.legend(); figure.tight_layout()
    figure.savefig(path, dpi=150); plt.close(figure)


def build_report(output, run_manifest, software_versions, classification, intro_images,
                 affinity_images, contour, individual, gmm, gmm_plots, adaptive):
    """Write HTML and its enabled report-local figures, without scientific work."""
    output = Path(output); output.parent.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(Path(run_manifest).read_text())
    software = _read_tsv(software_versions)
    sections = []
    overview = {key: manifest.get(key, '') for key in ('run_time', 'genome_build', 'populations', 'chromosomes', 'archaic_references', 'enabled_modules', 'git_commit')}
    sections.append('<h2>A. Run overview</h2>' + _table([overview]) + '<h3>Software versions</h3>' + _table(software))
    class_rows = _read_tsv(classification)
    class_plot = output.parent / 'classification_summary.png'
    _barplot(class_rows, 'genome_coverage_mb', class_plot, 'Classification genome coverage (Mb)')
    sections.append('<h2>B. Classification summary</h2>' + _table(class_rows) + '<img src="classification_summary.png" width="700">')
    sections.append('<h2>C. Introgression landscape</h2>' + (''.join(_image(path, 900) for path in intro_images) if intro_images else '<p>Not enabled.</p>'))
    sections.append('<h2>D. Landscape + affinity</h2>' + (''.join(_image(path, 900) for path in affinity_images) if affinity_images else '<p>Not enabled.</p>'))
    sections.append('<h2>E. Neanderthal × Denisovan contour plot</h2>' + (_image(contour, 700) if contour else '<p>Not enabled.</p>'))
    individual_rows = _read_tsv(individual)
    if individual_rows:
        individual_plot = output.parent / 'individual_summary.png'
        _barplot(individual_rows, 'per_individual_introgressed_mb', individual_plot, 'Per-individual archaic haplotype amount (Mb)')
        individual_html = _table(individual_rows) + '<img src="individual_summary.png" width="700">'
    else:
        individual_html = '<p>Not enabled.</p>'
    sections.append('<h2>F. Individual summary</h2>' + individual_html)
    gmm_rows = _read_tsv(gmm)
    selected = {f"{row['population']}__{row['reference']}" for row in gmm_rows
                if row.get('selected_components') not in ('NA', '', None) and int(row['selected_components']) > 1}
    figures = [path for path in gmm_plots if Path(path).stem.replace('.', '__') in selected]
    sections.append('<h2>G. GMM summary</h2>' + (_table(gmm_rows) + ''.join(_image(path, 500) for path in figures) if gmm_rows else '<p>Not enabled.</p>'))
    sections.append('<h2>H. Adaptive candidate summary</h2>' + (_table(_read_tsv(adaptive)) if adaptive else '<p>Not enabled.</p>'))
    page = '<!doctype html><html><head><meta charset="utf-8"><style>body{font-family:sans-serif}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px}img{display:block;margin:10px}</style></head><body><h1>SPrime workflow report</h1>' + ''.join(sections) + '</body></html>'
    output.write_text(page)


if 'snakemake' in globals():
    build_report(snakemake.output[0], snakemake.input.run, snakemake.input.software,
                 snakemake.input.classification, list(snakemake.input.intro),
                 list(snakemake.input.affinity), snakemake.input.contour or None,
                 snakemake.input.individual or None, snakemake.input.gmm or None,
                 list(snakemake.input.gmm_plots), snakemake.input.adaptive or None)
