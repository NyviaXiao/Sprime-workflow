"""Shared table I/O and dependency-free scientific analysis helpers.

Coordinates in output tables are 1-based and inclusive. BED masks are handled
by map_arch and remain 0-based, half-open. Missing match rates are serialized
as the literal ``NA`` so R and pandas can read the same files consistently.
"""
import csv
import gzip
from collections import defaultdict
from pathlib import Path

# Canonical schemas keep column order stable across populations and reruns.
BASE = ['population', 'chromosome', 'segment_id', 'start', 'end',
        'length_bp', 'sprime_score', 'marker_count']
CALLS = ['population', 'sample_id', 'haplotype', 'chromosome', 'start',
         'end', 'length_bp', 'marker_count', 'segment_id', 'archaic_class']


def op(path, mode='rt'):
    """Open plain text or gzip text using the filename suffix."""
    return gzip.open(path, mode, newline='') if str(path).endswith('.gz') else open(path, mode, newline='')


def read(path):
    """Read a headered TSV into dictionaries; intended for summary tables."""
    with op(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, fields, rows):
    """Write dictionaries as a deterministic TSV, optionally gzip-compressed."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with op(path, 'wt') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def collect(paths, output):
    """Concatenate header-compatible TSV files without loading all rows."""
    with op(paths[0]) as f:
        fields = next(csv.reader(f, delimiter='\t'))
    def rows():
        for path in paths:
            with op(path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                if reader.fieldnames != fields:
                    raise ValueError(f'Incompatible table columns: {path}')
                yield from reader
    write(output, fields, rows())


def score_rows(path):
    """Yield SPrime/mscore records using names rather than column positions."""
    # Historical files use a mixture of tabs and spaces, so split on arbitrary
    # whitespace while validating that every row matches the header width.
    with op(path) as f:
        head = f.readline().strip().split()
        required = {'CHROM', 'POS', 'SEGMENT', 'SCORE', 'ALLELE', 'REF', 'ALT'}
        if not required <= set(head):
            raise ValueError(f'{path}: missing score fields {required - set(head)}')
        for line in f:
            vals = line.split()
            if not vals:
                continue
            if len(vals) != len(head):
                raise ValueError(f'Malformed score row: {path}')
            yield dict(zip(head, vals))


def summarize(paths, pop, refs, minimum, score_gt=None):
    """Create one segment row with counts/rates for every ancient reference.

    ``score_gt`` reproduces pre_data.r's strict SCORE filter when preparing
    GMM input. It is deliberately independent of physical segment length.
    Match rate is match / (match + mismatch); notcomp is excluded.
    """
    segments = {}
    counts = defaultdict(lambda: [0, 0, 0])
    # For a chromosome, every independently annotated ancient reference must
    # describe the exact same underlying SPrime markers.
    reference_keys, previous_chrom = {}, None
    for path in sorted(paths, key=lambda p: (Path(p).name, str(p))):
        if Path(path).name != previous_chrom:
            reference_keys.clear()
            previous_chrom = Path(path).name
        rid = Path(path).parent.name
        tag = refs[rid]['tag']
        seen = set()
        for row in score_rows(path):
            key = (row['CHROM'], row['SEGMENT'])
            # Include alleles in the identity key to catch duplicated or
            # inconsistent records at the same genomic position.
            site = key + (row['POS'], row['REF'], row['ALT'], row['ALLELE'])
            if site in seen:
                raise ValueError(f'Duplicate score allele: {site}')
            seen.add(site)
            if score_gt is not None and float(row['SCORE']) <= score_gt:
                continue
            pos, score = int(row['POS']), float(row['SCORE'])
            # A segment's displayed boundaries are its first and last retained
            # SPrime markers, not inferred recombination breakpoints.
            if key not in segments:
                segments[key] = dict(zip(BASE, [pop, *key, pos, pos, 1, score, 0]))
            item = segments[key]
            item['start'], item['end'] = min(item['start'], pos), max(item['end'], pos)
            item['sprime_score'] = max(item['sprime_score'], score)
            state = row[tag]
            if state not in ('match', 'mismatch', 'notcomp'):
                raise ValueError(f'Unknown match state {state}')
            counts[key, rid][('match', 'mismatch', 'notcomp').index(state)] += 1
        # Each reference must annotate the identical underlying score rows.
        chrom = Path(path).name
        if chrom in reference_keys and seen != reference_keys[chrom]:
            raise ValueError(f'Mismatched reference score rows: {path}')
        reference_keys[chrom] = seen
    def order(k):
        """Sort numeric chromosomes naturally and named contigs afterwards."""
        chrom = k[0].removeprefix('chr')
        return (int(chrom) if chrom.isdigit() else 1000, k[0], segments[k]['start'], k[1])
    out = []
    for key in sorted(segments, key=order):
        item = segments[key]
        item['length_bp'] = item['end'] - item['start'] + 1
        for rid in refs:
            matched, mismatch, notcomp = counts[key, rid]
            # ``notcomp`` represents sites outside the callable ancient mask or
            # without an ancient genotype and therefore cannot enter the rate.
            callable_n = matched + mismatch
            item.update({f'{rid}_matched': matched, f'{rid}_mismatch': mismatch,
                         f'{rid}_callable': callable_n, f'{rid}_notcomp': notcomp,
                         f'{rid}_match_rate': matched / callable_n if callable_n >= minimum and callable_n else 'NA'})
            item['marker_count'] = max(item['marker_count'], callable_n + notcomp)
        out.append(item)
    return out


def summary_fields(refs):
    """Return the stable wide-table schema for the configured references."""
    return BASE + [f'{r}_{s}' for r in refs for s in ('matched', 'mismatch', 'callable', 'notcomp', 'match_rate')]


def rate(row, rid):
    """Parse a match rate, mapping the on-disk NA representation to None."""
    v = row[f'{rid}_match_rate']
    return None if v in ('NA', '', None) else float(v)


def classify(row, cfg):
    """Assign a class using the maximum available rate in each archaic group."""
    nean_rates = [rate(row, rid) for rid in cfg['neanderthal_references']]
    deni_rates = [rate(row, rid) for rid in cfg['denisovan_references']]
    nean_rates = [v for v in nean_rates if v is not None]
    deni_rates = [v for v in deni_rates if v is not None]
    if not nean_rates or not deni_rates:
        return 'ambiguous'
    max_nean, max_deni = max(nean_rates), max(deni_rates)
    is_nean = (max_nean > cfg['neanderthal_gt'] and
               max_deni < cfg['denisovan_lt_for_neanderthal'])
    is_deni = (max_deni > cfg['denisovan_gt'] and
               max_nean < cfg['neanderthal_lt_for_denisovan'])
    if is_nean:
        return 'neanderthal'
    if is_deni:
        return 'denisovan'
    return 'ambiguous'


def gmm_pass(row, cfg):
    """Apply group-level GMM eligibility using only callable references."""
    nr, dr = cfg['neanderthal_references'], cfg['denisovan_references']
    nean = [rate(row, r) for r in nr
            if int(row[f'{r}_callable']) >= cfg['min_callable'] and rate(row, r) is not None]
    deni = [rate(row, r) for r in dr
            if int(row[f'{r}_callable']) >= cfg['min_callable'] and rate(row, r) is not None]
    if not nean or not deni:
        return False
    return max(nean) < cfg['neanderthal_rate_lt'] and max(deni) > cfg['denisovan_rate_gt']


def gmm_target_pass(row, cfg, target):
    """Require one target reference to be callable for its own GMM input."""
    if int(row[f'{target}_callable']) < cfg['min_callable'] or rate(row, target) is None:
        return False
    return gmm_pass(row, cfg)


def classification_summary_rows(class_tables, wide_tables, populations):
    """Return one count/coverage row for each population and class."""
    output = []
    for population in populations:
        lengths = {(row['chromosome'], row['segment_id']): float(row['length_bp'])
                   for row in wide_tables[population]}
        for label in ('neanderthal', 'denisovan', 'ambiguous'):
            selected = [row for row in class_tables[population]
                        if row['archaic_class'] == label]
            output.append({'population': population, 'archaic_class': label,
                           'n_segments': len(selected),
                           'genome_coverage_mb':
                           sum(lengths[(row['chromosome'], row['segment_id'])]
                               for row in selected) / 1_000_000})
    return output


def individual_summary_rows(call_tables, sample_rows, populations):
    """Summarize haplotype tract lengths and unique target individuals."""
    target_counts = {population: sum(r['role'] == 'target' and r['population'] == population
                                    for r in sample_rows)
                     for population in populations}
    output = []
    for population in populations:
        rows = call_tables[population]
        for label in ('neanderthal', 'denisovan', 'ambiguous'):
            selected = [r for r in rows if r['archaic_class'] == label]
            output.append({'population': population, 'archaic_class': label,
                           'per_individual_introgressed_mb':
                           sum(float(r['length_bp']) for r in selected) /
                           target_counts[population] / 1_000_000,
                           'n_introgressed_individuals': len({r['sample_id'] for r in selected})})
    return output


def runs(markers, minimum=2):
    """Yield consecutive carried-marker runs.

    Consecutive means adjacent SPrime markers, not adjacent bases. False or a
    missing genotype terminates a run. Returned endpoints are 1-based marker
    coordinates and are not estimates of the biological tract breakpoints.
    """
    start = end = None
    count = 0
    for pos, carries in markers:
        if carries:
            if start is None:
                start = pos
            end = pos
            count += 1
        else:
            if count >= minimum:
                yield start, end, count
            start = end = None
            count = 0
    if count >= minimum:
        yield start, end, count
