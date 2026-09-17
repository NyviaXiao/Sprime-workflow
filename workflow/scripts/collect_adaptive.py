#!/usr/bin/env python3
"""Collect chromosome adaptive outputs and select group-specific Top2 candidates.

Chromosome-level adaptive screening is deliberately completed upstream.  This
module only combines those results, ranks Neanderthal- and Denisovan-associated
segments independently, and calculates within-group cross-population overlaps.
"""
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def rank_top2(rows, adaptive_group):
    """Return the two highest-core-AF rows passing one archaic group.

    A dual-pass segment is intentionally eligible for both rankings.  Copying
    its row prevents the group-specific rank from overwriting the other copy.
    """
    pass_field = f'{adaptive_group}_pass'
    ranked = sorted((dict(row) for row in rows if int(row[pass_field]) == 1),
                    key=lambda row: (-float(row['core_mean_AF']), str(row['chromosome']),
                                     int(row['candidate_start']), row['segment_id']))
    for rank, row in enumerate(ranked, 1):
        row['adaptive_group'] = adaptive_group
        row['adaptive_rank'] = rank
        row['selected_top2'] = int(rank <= 2)
    return ranked[:2]


def shared_intervals(rows):
    """Return BED spans shared by two populations, never across archaic groups."""
    intervals = defaultdict(list)
    for row in rows:
        key = (row['adaptive_group'], str(row['chromosome']))
        intervals[key].append((int(row['candidate_start']), int(row['candidate_end']), row['population']))
    result = []
    for (adaptive_group, chromosome), items in intervals.items():
        points = sorted({point for start, end, _ in items for point in (start, end)})
        for start, end in zip(points, points[1:]):
            populations = sorted({pop for left, right, pop in items if left < end and right > start})
            if len(populations) < 2:
                continue
            names = ','.join(populations)
            if (result and result[-1]['adaptive_group'] == adaptive_group and
                    result[-1]['chromosome'] == chromosome and result[-1]['end'] == start and
                    result[-1]['populations'] == names):
                result[-1]['end'] = end
            else:
                result.append({'adaptive_group': adaptive_group, 'chromosome': chromosome,
                               'start': start, 'end': end,
                               'populations': names, 'n_populations': len(populations)})
    return result


def _read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def _write(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def collect(inputs, outputs):
    variants, core_rows, segments = [], [], []
    for path in inputs:
        if str(path).endswith('.segments.tsv'):
            segments.extend(_read(path))
        elif str(path).endswith('.core.tsv.gz'):
            core_rows.extend(_read(path))
        else:
            variants.extend(_read(path))
    for row in segments:
        row['pass_flag'] = int(row['pass_flag'])
    selected = []
    for population in sorted({row['population'] for row in segments}):
        population_rows = [row for row in segments if row['population'] == population]
        for adaptive_group in ('neanderthal', 'denisovan'):
            selected.extend(rank_top2(population_rows, adaptive_group))
    fields = list(segments[0]) if segments else ['population','chromosome','segment_id','candidate_start','candidate_end','max_AF','core_mean_AF','core_variant_count','neanderthal_pass','denisovan_pass','pass_flag']
    for name in ('adaptive_group', 'adaptive_rank', 'selected_top2'):
        if name not in fields:
            fields.append(name)
    _write(outputs['variants'], list(variants[0]) if variants else ['population'], variants)
    _write(outputs['core'], list(core_rows[0]) if core_rows else ['population'], core_rows)
    _write(outputs['segments'], fields, segments)
    _write(outputs['top2'], fields, selected)
    _write(outputs['shared'], ['adaptive_group','chromosome','start','end','populations','n_populations'], shared_intervals(selected))


if 'snakemake' in globals():
    collect(snakemake.input, dict(snakemake.output))
