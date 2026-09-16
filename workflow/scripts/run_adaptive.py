#!/usr/bin/env python3
"""Adaptive-introgression filtering on SPrime markers.

The small pure functions below are deliberately independent of Snakemake so
the AF, pruning, threshold, ranking and interval semantics can be tested
without a VCF fixture.  The script mode joins those functions to one
population/chromosome and a final collection job.
"""
import csv
import gzip
import hashlib
from collections import defaultdict
from pathlib import Path


def introgressed_af(alt_af, sprime_allele):
    return float(alt_af) if int(sprime_allele) == 1 else 1.0 - float(alt_af)


def core_variants(rows, min_af=0.30, max_drop=0.20):
    eligible = [r for r in rows if float(r['introgressed_AF']) > min_af]
    if not eligible:
        return [], None
    maximum = max(float(r['introgressed_AF']) for r in eligible)
    return [r for r in eligible if float(r['introgressed_AF']) >= maximum - max_drop], maximum


def adaptive_pass(callable_sites, match_rate, group, min_callable=10,
                  nean_threshold=0.50, deni_threshold=0.40):
    if int(callable_sites) < min_callable or match_rate in (None, 'NA'):
        return False
    threshold = nean_threshold if group == 'neanderthal' else deni_threshold
    return float(match_rate) >= threshold


def rank_top2(rows, top_regions=2):
    ranked = sorted((r for r in rows if r['pass_flag']),
                    key=lambda r: (-float(r['core_mean_AF']), str(r['chromosome']),
                                   int(r['candidate_start']), r['segment_id']))
    for i, row in enumerate(ranked, 1):
        row['adaptive_rank'] = i
        row['selected_top2'] = '1' if i <= top_regions else '0'
    return ranked[:top_regions]


def shared_intervals(rows):
    """Return all spans covered by at least two populations (half-open BED)."""
    by_chrom = defaultdict(list)
    for row in rows:
        by_chrom[str(row['chromosome'])].append((int(row['start']), int(row['end']), row['population']))
    out = []
    for chrom, intervals in by_chrom.items():
        points = sorted({p for s, e, _ in intervals for p in (s, e)})
        for left, right in zip(points, points[1:]):
            active = sorted({pop for s, e, pop in intervals if s < right and e > left})
            if len(active) >= 2:
                if out and out[-1]['chromosome'] == chrom and out[-1]['end'] == left and out[-1]['populations'] == ','.join(active):
                    out[-1]['end'] = right
                else:
                    out.append({'chromosome': chrom, 'start': left, 'end': right,
                                'populations': ','.join(active), 'n_populations': len(active)})
    return out


def _read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        yield from csv.DictReader(f, delimiter='\t')

def _read_score(path):
    """Read SPrime's whitespace-delimited score format by header name."""
    with open(path) as f:
        head = f.readline().split()
        for line in f:
            values = line.split()
            if values:
                if len(values) != len(head): raise ValueError(f'Malformed score row: {path}')
                yield dict(zip(head, values))


def _open_write(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)


def collect(chrom_files, output_dir, top_regions):
    variants, cores, segments = [], [], []
    for path in chrom_files:
        if str(path).endswith('.segments.tsv'):
            segments.extend(_read(path))
        elif str(path).endswith('.core.tsv.gz'):
            cores.extend(_read(path))
        else:
            variants.extend(_read(path))
    for row in segments:
        row['pass_flag'] = str(row.get('pass_flag', '0'))
    top = []
    for pop in sorted({r['population'] for r in segments}):
        top.extend(rank_top2([r for r in segments if r['population'] == pop], top_regions))
    seg_fields = list(segments[0]) if segments else ['population']
    for f in ('adaptive_rank', 'selected_top2'):
        if f not in seg_fields: seg_fields.append(f)
    _open_write(Path(output_dir)/'adaptive/variants.tsv.gz', list(variants[0]) if variants else ['population'], variants)
    _open_write(Path(output_dir)/'adaptive/core_variants.tsv.gz', list(cores[0]) if cores else ['population'], cores)
    _open_write(Path(output_dir)/'adaptive/segment_summary.tsv', seg_fields, segments)
    _open_write(Path(output_dir)/'adaptive/top2_candidates.tsv', seg_fields, top)
    _open_write(Path(output_dir)/'adaptive/shared_top2_regions.tsv', ['chromosome','start','end','populations','n_populations'], shared_intervals(top))


if 'snakemake' in globals():
    mode = snakemake.params.mode
    if mode == 'collect':
        collect(list(snakemake.input), snakemake.params.outdir, int(snakemake.params.top_regions))
    else:
        # Marker-level extraction is intentionally explicit: upstream SPrime
        # allele and archaic status are never replaced by ordinary ALT AF.
        import pysam
        pop, chrom = snakemake.wildcards.population, snakemake.wildcards.chrom
        cfg = snakemake.params.adaptive
        score_rows = list(_read_score(snakemake.input.score))
        samples = [r['sample_id'] for r in _read(snakemake.input.samples) if r['role'] == 'target' and r['population'] == pop]
        by_key = {}
        with pysam.VariantFile(snakemake.input.vcf) as vf:
            for rec in vf.fetch(chrom):
                if not rec.alts or len(rec.alts) != 1: continue
                vals = [allele for sample in samples
                        for allele in (rec.samples[sample].get('GT') or ())
                        if allele is not None]
                if vals: by_key[(rec.pos, rec.ref, rec.alts[0])] = sum(a == 1 for a in vals) / len(vals)
        status = {}
        for path, ref in zip(snakemake.input.mscore, snakemake.params.references):
            for r in _read(path): status[(int(r['POS']), r['REF'], r['ALT'], ref)] = r.get(snakemake.params.tags[ref], 'notcomp')
        marker_rows = []
        for r in score_rows:
            key = (int(r['POS']), r['REF'], r['ALT'])
            if key not in by_key: continue
            af = introgressed_af(by_key[key], r['ALLELE'])
            if af > float(cfg['min_allele_frequency']): marker_rows.append({'population':pop,'chromosome':chrom,'position':r['POS'],'segment_id':r['SEGMENT'],'sprime_allele':r['ALLELE'],'introgressed_AF':af, **{f'{ref}_status': status.get((*key, ref), 'notcomp') for ref in snakemake.params.references}})
        by_seg = defaultdict(list)
        for r in marker_rows: by_seg[r['segment_id']].append(r)
        seg_rows = []
        for seg, rows in by_seg.items():
            core, maximum = core_variants(rows, float(cfg['min_allele_frequency']), float(cfg['max_frequency_drop']))
            if not core: continue
            out = {'population':pop,'chromosome':chrom,'segment_id':seg,'candidate_start':min(int(r['position'])-1 for r in core),'candidate_end':max(int(r['position']) for r in core),'max_AF':maximum,'core_mean_AF':sum(float(r['introgressed_AF']) for r in core)/len(core),'core_variant_count':len(core)}
            for ref in snakemake.params.references:
                st = [r[f'{ref}_status'] for r in core]; call = sum(x in ('match','mismatch') for x in st); mat = sum(x == 'match' for x in st); rate = mat/call if call else 'NA'; out[f'{ref}_adaptive_callable']=call; out[f'{ref}_adaptive_match_rate']=rate; out[f'{ref}_adaptive_pass']=int(adaptive_pass(call, rate, snakemake.params.groups[ref], int(cfg['min_callable_sites']), float(cfg['neanderthal_match_threshold']), float(cfg['denisovan_match_threshold'])))
            out['neanderthal_pass']=int(any(out[f'{r}_adaptive_pass'] for r in snakemake.params.neanderthal_refs)); out['denisovan_pass']=int(any(out[f'{r}_adaptive_pass'] for r in snakemake.params.denisovan_refs)); out['pass_flag']=int(out['neanderthal_pass'] or out['denisovan_pass']); seg_rows.append(out); cores.extend(core)
        base = Path(snakemake.output.segments); _open_write(base, list(seg_rows[0]) if seg_rows else ['population'], seg_rows); _open_write(snakemake.output.variants, list(marker_rows[0]) if marker_rows else ['population'], marker_rows); _open_write(snakemake.output.core, list(cores[0]) if cores else ['population'], cores)
