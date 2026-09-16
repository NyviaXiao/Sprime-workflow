#!/usr/bin/env python3
"""Build a compact HTML presentation of final artifacts."""
import csv, html
from pathlib import Path

SECTION_ORDER = ['Run overview','Classification summary','Introgression landscape','Landscape + affinity','Neanderthal × Denisovan contour plot','Individual summary','GMM summary','Adaptive candidate summary']

def _read(path):
    with open(path, newline='') as f: return list(csv.DictReader(f, delimiter='\t'))

def _table(rows):
    if not rows: return '<p>No rows.</p>'
    fields=list(rows[0]); out=['<table><thead><tr>'+''.join(f'<th>{html.escape(x)}</th>' for x in fields)+'</tr></thead><tbody>']
    for r in rows: out.append('<tr>'+''.join(f'<td>{html.escape(str(r.get(x,"")))}</td>' for x in fields)+'</tr>')
    return ''.join(out)+'</tbody></table>'

def _barplot(rows, value, path, title):
    """Create the two compact grouped summaries required by the report."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    pops=sorted({r['population'] for r in rows}); classes=['neanderthal','denisovan','ambiguous']
    fig, ax=plt.subplots(figsize=(8,4)); width=.24
    for i, cl in enumerate(classes):
        vals=[float(next((r[value] for r in rows if r['population']==p and r['archaic_class']==cl), 0)) for p in pops]
        ax.bar([x+(i-1)*width for x in range(len(pops))], vals, width, label=cl)
    ax.set_xticks(range(len(pops))); ax.set_xticklabels(pops); ax.set_ylabel(value); ax.set_title(title); ax.legend(); fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)

def build(output, run, classification, individual, adaptive, images, contour, gmm, gmm_plots=None):
    sections=[]
    cfg=_read(run) if str(run).endswith('.tsv') else None
    sections.append(('<h2>A. Run overview</h2><p>Run metadata and software versions are stored in provenance/. Full resolved configuration is not embedded here.</p>'))
    class_rows=_read(classification); _barplot(class_rows, 'genome_coverage_mb', Path(output).parent/'classification_summary.png', 'Classification genome coverage (Mb)')
    sections.append('<h2>B. Classification summary</h2>'+_table(class_rows)+'<img src="classification_summary.png" width="700">')
    sections.append('<h2>C. Introgression landscape</h2>'+''.join(f'<img src="../{html.escape(str(p))}" width="900">' for p in images.get('intro',[])))
    sections.append('<h2>D. Landscape + affinity</h2>'+''.join(f'<img src="../{html.escape(str(p))}" width="900">' for p in images.get('affinity',[])))
    sections.append('<h2>E. Neanderthal × Denisovan contour plot</h2>'+ (f'<img src="../{html.escape(str(contour))}" width="700">' if contour else '<p>Not enabled.</p>'))
    ind_rows=_read(individual) if individual else []; ind_plot=''
    if ind_rows:
        _barplot(ind_rows, 'per_individual_introgressed_mb', Path(output).parent/'individual_summary.png', 'Per-individual archaic haplotype amount (Mb)'); ind_plot='<img src="individual_summary.png" width="700">'
    sections.append('<h2>F. Individual summary</h2>'+(_table(ind_rows)+ind_plot if ind_rows else '<p>Not enabled.</p>'))
    g_rows=_read(gmm) if gmm else []; chosen={r['population']+'__'+r['reference'] for r in g_rows if r.get('selected_components','0') not in ('NA','0','1') and int(r['selected_components']) > 1}
    figs=''.join(f'<img src="../{html.escape(str(p))}" width="500">' for p in (gmm_plots or []) if Path(p).stem.replace('.','__') in chosen)
    sections.append('<h2>G. GMM summary</h2>'+(_table(g_rows)+figs if g_rows else '<p>Not enabled.</p>'))
    sections.append('<h2>H. Adaptive candidate summary</h2>'+(_table(_read(adaptive)) if adaptive else '<p>Not enabled.</p>'))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text('<!doctype html><html><head><meta charset="utf-8"><style>body{font-family:sans-serif}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:4px}img{display:block;margin:10px}</style></head><body><h1>SPrime workflow report</h1>'+''.join(sections)+'</body></html>')

if 'snakemake' in globals():
    build(snakemake.output.html, snakemake.input.run, snakemake.input.classification, snakemake.input.individual, snakemake.input.adaptive if snakemake.params.adaptive else None, {'intro':snakemake.input.intro,'affinity':snakemake.input.affinity}, snakemake.input.contour, snakemake.input.gmm, snakemake.input.gmm_plots)
