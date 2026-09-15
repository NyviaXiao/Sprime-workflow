"""Validate inputs before expensive jobs and record run provenance.

The checks focus on assumptions made by SPrime and the bundled map_arch:
sample membership, chromosome naming, one ancient individual per VCF and
sorted single-chromosome BED masks.
"""
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from core import read, write, op


def validate(s, c):
    """Run fail-fast validation and write QC/provenance output files."""
    import pysam
    samples = read(c['samples'])
    wanted = {r['sample_id'] for r in samples}
    records, modern = [], {}
    # Validate modern and ancient inputs chromosome by chromosome. This scans
    # ancient records intentionally because map_arch indexes solely by POS.
    for chrom in c['chromosomes']:
        path = c['vcf'].format(chrom=chrom)
        with pysam.VariantFile(path) as v:
            present = set(v.header.samples)
            if not wanted <= present:
                raise ValueError(f'{path}: missing samples {wanted-present}')
            if chrom not in v.header.contigs:
                raise ValueError(f'{path}: chromosome {chrom} missing in header')
            modern[path] = True
        records.append({'check': 'modern_header', 'detail': f'{chrom}: {path}', 'status': 'PASS'})
        for ref in c['archaic_references']:
            av = ref['vcf'].format(chrom=chrom)
            with pysam.VariantFile(av) as v:
                if len(v.header.samples) != 1:
                    raise ValueError(f'{av}: map_arch requires exactly one ancient sample')
                # Bundled map_arch indexes by POS, so multi-chromosome inputs are unsafe.
                for record in v:
                    if record.contig != chrom:
                        raise ValueError(f'{av}: must contain only chromosome {chrom}')
            # map_arch treats BED intervals as 0-based, half-open callable
            # regions. Sorted, non-overlapping intervals avoid ambiguous masks.
            mask = ref['mask'].format(chrom=chrom)
            last_end = -1
            with op(mask) as f:
                for line in f:
                    if not line.strip() or line.startswith(('#', 'track', 'browser')):
                        continue
                    fields = line.split()
                    if len(fields) < 3 or fields[0] != chrom:
                        raise ValueError(f'{mask}: must contain only chromosome {chrom}')
                    start, end = map(int, fields[1:3])
                    if start < 0 or end <= start or start < last_end:
                        raise ValueError(f'{mask}: BED must be sorted, nonoverlapping, 0-based half-open')
                    last_end = end
            records.append({'check': 'archaic', 'detail': f'{ref["id"]}: {chrom}', 'status': 'PASS'})
    # SPrime expects chromosome identifiers in its genetic map to match those
    # supplied through chrom= exactly (for example, 1 versus chr1).
    map_chroms = set()
    with open(c['genetic_map']) as f:
        for line in f:
            fields = line.split()
            if fields:
                map_chroms.add(fields[0])
    if not set(c['chromosomes']) <= map_chroms:
        raise ValueError('Genetic map chromosome names do not match configuration')
    for pop in c['populations']:
        count = sum(r['role'] == 'target' and r['population'] == pop for r in samples)
        records.append({'check': 'target_samples', 'detail': f'{pop}: {count}', 'status': 'PASS'})
    records.append({'check': 'genome_build', 'detail': c['genome_build'] + ': user-declared, not inferred', 'status': 'DECLARED'})
    # Capture complete version strings because tools often include build and
    # linked-library details after the first line.
    versions = {}
    for tool, args in [('bcftools', ['bcftools', '--version']), ('java', ['java', '-version']), ('R', ['Rscript', '--version'])]:
        r = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        versions[tool] = r.stdout.strip()
    def digest(path):
        """SHA-256 small, essential inputs without loading them into memory."""
        h = hashlib.sha256()
        with open(path, 'rb') as f:
            for block in iter(lambda: f.read(1048576), b''):
                h.update(block)
        return h.hexdigest()
    # Large VCFs are represented by path, size and mtime to keep startup cost
    # reasonable; small code/tool inputs receive a cryptographic checksum.
    info = {'time_utc': datetime.now(timezone.utc).isoformat(),
            'resolved_config': str(Path(c['outdir']) / 'resolved_config.json'), 'versions': versions,
            'sha256': {k: digest(c[k]) for k in ('samples', 'sprime_jar', 'map_arch')},
            'inputs': [{'path': str(p), 'size': Path(p).stat().st_size, 'mtime_ns': Path(p).stat().st_mtime_ns} for p in s.input]}
    Path(s.output[1]).write_text(json.dumps(info, indent=2))
    write(s.output[0], ['check', 'detail', 'status'], records)
