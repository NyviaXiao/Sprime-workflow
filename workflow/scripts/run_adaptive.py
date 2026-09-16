#!/usr/bin/env python3
"""Per-chromosome adaptive filtering from SPrime markers and mscore states."""
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def introgressed_af(alt_af, sprime_allele):
    """Orient ALT frequency to SPrime's inferred allele (0=REF, 1=ALT)."""
    return float(alt_af) if int(sprime_allele) == 1 else 1.0 - float(alt_af)


def core_variants(rows, min_af, max_drop):
    """Apply strict AF filtering, then retain sites within max_drop of max AF."""
    eligible = [row for row in rows if float(row['introgressed_AF']) > min_af]
    if not eligible:
        return [], None
    maximum = max(float(row['introgressed_AF']) for row in eligible)
    return [row for row in eligible if float(row['introgressed_AF']) >= maximum - max_drop], maximum


def adaptive_pass(callable_sites, match_rate, group, min_callable, nean_threshold, deni_threshold):
    if int(callable_sites) < min_callable or match_rate == 'NA':
        return False
    threshold = nean_threshold if group == 'neanderthal' else deni_threshold
    return float(match_rate) >= threshold


def build_chromosome_rows(population, chromosome, score_rows, alt_af, states, cfg,
                          references, groups, neanderthal_refs, denisovan_refs):
    """Return marker/core/segment rows for one population and chromosome."""
    markers = []
    for score in score_rows:
        key = (int(score['POS']), score['REF'], score['ALT'])
        if key not in alt_af:
            continue
        af = introgressed_af(alt_af[key], score['ALLELE'])
        row = {'population': population, 'chromosome': chromosome,
               'position': score['POS'], 'segment_id': score['SEGMENT'],
               'sprime_allele': score['ALLELE'], 'introgressed_AF': af}
        row.update({f'{ref}_status': states[ref].get(key, 'notcomp') for ref in references})
        markers.append(row)
    by_segment = defaultdict(list)
    for row in markers:
        by_segment[row['segment_id']].append(row)
    core_rows, segments = [], []
    for segment_id, rows in by_segment.items():
        core, maximum = core_variants(rows, float(cfg['min_allele_frequency']), float(cfg['max_frequency_drop']))
        if not core:
            continue
        core_rows.extend(core)
        segment = {'population': population, 'chromosome': chromosome, 'segment_id': segment_id,
                   'candidate_start': min(int(row['position']) - 1 for row in core),
                   'candidate_end': max(int(row['position']) for row in core),
                   'max_AF': maximum,
                   'core_mean_AF': sum(float(row['introgressed_AF']) for row in core) / len(core),
                   'core_variant_count': len(core)}
        for ref in references:
            ref_states = [row[f'{ref}_status'] for row in core]
            callable_sites = sum(state in ('match', 'mismatch') for state in ref_states)
            matches = sum(state == 'match' for state in ref_states)
            match_rate = matches / callable_sites if callable_sites else 'NA'
            segment[f'{ref}_adaptive_callable'] = callable_sites
            segment[f'{ref}_adaptive_match_rate'] = match_rate
            segment[f'{ref}_adaptive_pass'] = int(adaptive_pass(
                callable_sites, match_rate, groups[ref], int(cfg['min_callable_sites']),
                float(cfg['neanderthal_match_threshold']), float(cfg['denisovan_match_threshold'])))
        segment['neanderthal_pass'] = int(any(segment[f'{ref}_adaptive_pass'] for ref in neanderthal_refs))
        segment['denisovan_pass'] = int(any(segment[f'{ref}_adaptive_pass'] for ref in denisovan_refs))
        segment['pass_flag'] = int(segment['neanderthal_pass'] or segment['denisovan_pass'])
        segments.append(segment)
    return markers, core_rows, segments


def _read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def _read_score(path):
    with open(path) as handle:
        header = handle.readline().split()
        for line in handle:
            values = line.split()
            if values:
                if len(values) != len(header):
                    raise ValueError(f'Malformed score row: {path}')
                yield dict(zip(header, values))


def _write(path, fields, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)


def _alt_frequencies(vcf_path, chromosome, samples):
    import pysam
    values = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.fetch(chromosome):
            if not record.alts or len(record.alts) != 1:
                continue
            alleles = [allele for sample in samples for allele in (record.samples[sample].get('GT') or ()) if allele is not None]
            if alleles:
                values[(record.pos, record.ref, record.alts[0])] = sum(a == 1 for a in alleles) / len(alleles)
    return values


if 'snakemake' in globals():
    population, chromosome = snakemake.wildcards.population, snakemake.wildcards.chrom
    cfg = snakemake.params.adaptive
    samples = [row['sample_id'] for row in _read_table(snakemake.input.samples)
               if row['role'] == 'target' and row['population'] == population]
    states = {}
    for path, ref in zip(snakemake.input.mscore, snakemake.params.references):
        tag = snakemake.params.tags[ref]
        states[ref] = {(int(row['POS']), row['REF'], row['ALT']): row.get(tag, 'notcomp') for row in _read_table(path)}
    marker_rows, core_rows, segment_rows = build_chromosome_rows(
        population, chromosome, list(_read_score(snakemake.input.score)),
        _alt_frequencies(snakemake.input.vcf, chromosome, samples), states, cfg,
        snakemake.params.references, snakemake.params.groups,
        snakemake.params.neanderthal_refs, snakemake.params.denisovan_refs)
    marker_fields = list(marker_rows[0]) if marker_rows else ['population','chromosome','position','segment_id','sprime_allele','introgressed_AF']
    core_fields = list(core_rows[0]) if core_rows else marker_fields
    segment_fields = list(segment_rows[0]) if segment_rows else ['population','chromosome','segment_id','candidate_start','candidate_end','max_AF','core_mean_AF','core_variant_count','neanderthal_pass','denisovan_pass','pass_flag']
    _write(snakemake.output.variants, marker_fields, marker_rows)
    _write(snakemake.output.core, core_fields, core_rows)
    _write(snakemake.output.segments, segment_fields, segment_rows)
