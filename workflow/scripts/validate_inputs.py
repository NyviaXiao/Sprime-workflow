"""Lightweight fail-fast validation of essential workflow inputs."""

import json
from datetime import datetime, timezone
from pathlib import Path

from core import read, write


def validate(s, c):
    """Validate cheap input assumptions required by downstream tools."""
    import pysam

    samples = read(c['samples'])
    wanted_samples = {row['sample_id'] for row in samples}
    records = []

    # Modern VCF:
    # - required samples must exist
    # - chromosome naming must agree with config (e.g. 1 vs chr1)
    checked_modern = {}

    for chrom in c['chromosomes']:
        path = c['vcf'].format(chrom=chrom)

        if path not in checked_modern:
            with pysam.VariantFile(path) as vcf:
                checked_modern[path] = {
                    'samples': set(vcf.header.samples),
                    'contigs': set(vcf.header.contigs),
                }

        header = checked_modern[path]

        missing = wanted_samples - header['samples']
        if missing:
            raise ValueError(
                f'{path}: missing samples {sorted(missing)}'
            )

        if chrom not in header['contigs']:
            raise ValueError(
                f'{path}: chromosome {chrom} missing from VCF header'
            )

        records.append({
            'check': 'modern_vcf_header',
            'detail': f'{chrom}: {path}',
            'status': 'PASS',
        })

    # Ancient VCF:
    # - exactly one ancient sample
    # - chromosome naming must agree with config
    #
    # Header only: do NOT scan variant records.
    checked_ancient = {}

    for chrom in c['chromosomes']:
        for ref in c['archaic_references']:
            path = ref['vcf'].format(chrom=chrom)

            if path not in checked_ancient:
                with pysam.VariantFile(path) as vcf:
                    checked_ancient[path] = {
                        'n_samples': len(vcf.header.samples),
                        'contigs': set(vcf.header.contigs),
                    }

            header = checked_ancient[path]

            if header['n_samples'] != 1:
                raise ValueError(
                    f'{path}: map_arch requires exactly one ancient sample'
                )

            if chrom not in header['contigs']:
                raise ValueError(
                    f'{path}: chromosome {chrom} missing from VCF header'
                )

            records.append({
                'check': 'ancient_vcf_header',
                'detail': f'{ref["id"]}: {chrom}',
                'status': 'PASS',
            })

    # Genetic map chromosome naming must match config exactly.
    # This is a text file and is normally small compared with VCF inputs.
    map_chromosomes = set()

    with open(c['genetic_map']) as handle:
        for line in handle:
            fields = line.split()
            if fields:
                map_chromosomes.add(fields[0])

    missing_map_chromosomes = (
        set(c['chromosomes']) - map_chromosomes
    )

    if missing_map_chromosomes:
        raise ValueError(
            'Genetic map chromosome names do not match configuration: '
            f'missing {sorted(missing_map_chromosomes)}'
        )

    records.append({
        'check': 'genetic_map_chromosomes',
        'detail': ','.join(c['chromosomes']),
        'status': 'PASS',
    })

    records.append({
        'check': 'genome_build',
        'detail': (
            f'{c["genome_build"]}: user-declared, not inferred'
        ),
        'status': 'DECLARED',
    })

    # Minimal run information required by downstream provenance generation.
    info = {
        'time_utc': datetime.now(timezone.utc).isoformat(),
        'resolved_config': str(
            Path(c['outdir']) / 'resolved_config.json'
        ),
    }

    Path(s.output[1]).write_text(
        json.dumps(info, indent=2) + '\n'
    )

    write(
        s.output[0],
        ['check', 'detail', 'status'],
        records,
    )