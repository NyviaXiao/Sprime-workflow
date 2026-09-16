#!/usr/bin/env python3
"""Create one tidy affinity table for one population/reference pair.

The table is a view of the existing wide segment summary; no VCF or map_arch
work is repeated here.  Keeping this as a small pure-Python script also makes
the output useful to plotting and downstream reporting code.
"""
import csv
import gzip
from pathlib import Path


FIELDS = ['population', 'chromosome', 'segment_id', 'start', 'end', 'length_bp',
          'archaic_class', 'reference', 'match', 'mismatch', 'callable',
          'notcomp', 'match_rate']


def _rows(path):
    with open(path, newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def build(wide_path, classification_path, population, reference, output):
    classes = {(r['chromosome'], r['segment_id']): r['archaic_class']
               for r in _rows(classification_path)}
    rows = []
    prefix = f'{reference}_'
    for src in _rows(wide_path):
        key = (src['chromosome'], src['segment_id'])
        if key not in classes:
            raise ValueError(f'Missing classification for {key}')
        row = {'population': population, 'chromosome': src['chromosome'],
               'segment_id': src['segment_id'], 'start': src['start'],
               'end': src['end'], 'length_bp': src['length_bp'],
               'archaic_class': classes[key], 'reference': reference}
        for name in ('matched', 'mismatch', 'callable', 'notcomp', 'match_rate'):
            row[name if name != 'matched' else 'match'] = src[prefix + name]
        rows.append(row)
    def order(r):
        chrom = str(r['chromosome']).removeprefix('chr')
        return (int(chrom) if chrom.isdigit() else 1000, chrom, int(r['start']), r['segment_id'])
    rows.sort(key=order)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


if 'snakemake' in globals():
    build(snakemake.input.wide, snakemake.input.classification,
          snakemake.wildcards.population, snakemake.wildcards.ref,
          snakemake.output[0])
