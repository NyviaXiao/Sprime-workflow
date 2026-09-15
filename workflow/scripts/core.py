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
    """Assign mutually exclusive Nean/Deni labels; every remainder is ambiguous."""
    n, d = rate(row, cfg['neanderthal_reference']), rate(row, cfg['denisovan_reference'])
    if n is None or d is None:
        return 'ambiguous'
    nean = n > cfg['neanderthal_gt'] and d < cfg['denisovan_lt_for_neanderthal']
    deni = d > cfg['denisovan_gt'] and n < cfg['neanderthal_lt_for_denisovan']
    if nean:
        return 'neanderthal'
    if deni:
        return 'denisovan'
    return 'ambiguous'


def gmm_pass(row, cfg):
    """Apply the configurable pre_data-style allele and match-rate filters."""
    nr, dr = cfg['neanderthal_references'], cfg['denisovan_references']
    required = set(nr + dr + cfg['target_references'])
    if any(int(row[f'{r}_callable']) < cfg['min_callable'] or rate(row, r) is None for r in required):
        return False
    return (all(rate(row, r) < cfg['neanderthal_rate_lt'] for r in nr)
            and any(rate(row, r) > cfg['denisovan_rate_gt'] for r in dr))


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
