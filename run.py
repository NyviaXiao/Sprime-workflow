#!/usr/bin/env python3
"""Command-line entry point for the SPrime workflow.

This module validates the configuration, resolves relative paths and starts
Snakemake. Analysis code lives under ``workflow/``; this file deliberately
contains no scientific calculations.
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def resolve_config(path, overrides):
    """Load, normalize and validate configuration before Snakemake sees it.

    Paths written in YAML are relative to the YAML file. Paths supplied on the
    command line are relative to the caller's current working directory.
    The returned dictionary contains only absolute paths.
    """
    import yaml
    path = Path(path).resolve()
    c = yaml.safe_load(path.read_text())
    # Resolve top-level resources once so rules do not depend on their launch
    # directory. ``vcf`` may still contain the literal {chrom} placeholder.
    for key in ('vcf', 'samples', 'outdir', 'genetic_map', 'sprime_jar', 'map_arch'):
        value = overrides.get(key) or c[key]
        base = Path.cwd() if overrides.get(key) else path.parent
        c[key] = str((base / value).resolve())
    for ref in c['archaic_references']:
        for key in ('vcf', 'mask'):
            ref[key] = str((path.parent / ref[key]).resolve())
    # IDs become path components and column prefixes; tags become mscore
    # headers. Restrict both to portable, unambiguous identifiers.
    refs = c['archaic_references']
    ids = [r['id'] for r in refs]
    tags = [r['tag'] for r in refs]
    if not refs or len(set(ids)) != len(ids) or len(set(tags)) != len(tags):
        raise ValueError('Ancient reference IDs/tags must be nonempty and unique')
    for value in ids + tags:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,30}', value) or '__' in value:
            raise ValueError(f'Invalid reference ID/tag: {value}')
    if any(r['group'] not in ('neanderthal', 'denisovan') for r in refs):
        raise ValueError('group must be neanderthal or denisovan')
    c['chromosomes'] = [str(x) for x in c['chromosomes']]
    if not c['chromosomes'] or len(set(c['chromosomes'])) != len(c['chromosomes']):
        raise ValueError('chromosomes must be nonempty and unique')
    for chrom in c['chromosomes']:
        if not re.fullmatch(r'[A-Za-z0-9_]+', chrom):
            raise ValueError(f'Invalid chromosome: {chrom}')
    # Population names are derived from the sample sheet instead of being
    # duplicated in config.yaml.
    with open(c['samples'], newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if not rows or not {'sample_id', 'population', 'role'} <= rows[0].keys():
        raise ValueError('samples.tsv requires sample_id, population, role')
    if len({r['sample_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate sample IDs')
    if any(not r['sample_id'] or r['role'] not in ('target', 'outgroup') for r in rows):
        raise ValueError('Invalid sample ID/role')
    c['populations'] = sorted({r['population'] for r in rows if r['role'] == 'target'})
    if not c['populations'] or not any(r['role'] == 'outgroup' for r in rows):
        raise ValueError('At least one target population and outgroup sample required')
    if any(not re.fullmatch(r'[A-Za-z0-9_-]+', p) for p in c['populations']):
        raise ValueError('Population names may contain letters, digits, underscore, dash')
    # Every reference named by a downstream module must exist and must have
    # the expected biological group.
    cl = c['classification']
    for key in ('neanderthal_references', 'denisovan_references'):
        values = cl.get(key)
        if not isinstance(values, list) or not values:
            raise ValueError(f'classification.{key} must be a nonempty list')
        if len(set(values)) != len(values):
            raise ValueError(f'Duplicate reference in classification.{key}')
    used = cl['neanderthal_references'] + cl['denisovan_references']
    g = c['gmm']
    for name, value in [('summary_min_callable', c['summary_min_callable']),
                        ('gmm.min_callable', g['min_callable']), ('gmm.min_segments', g['min_segments']),
                        ('gmm.n_init', g['n_init']), ('individual.min_markers', c['individual']['min_markers']),
                        *c['resources'].items()]:
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise ValueError(f'{name} must be a positive integer')
    if g['enabled']:
        used += g['target_references'] + g['neanderthal_references'] + g['denisovan_references']
        if not all(g[k] for k in ('target_references', 'neanderthal_references', 'denisovan_references')):
            raise ValueError('GMM reference lists cannot be empty')
        if g['method'] not in ('legacy_chi_square', 'parametric_bootstrap'):
            raise ValueError('Unknown GMM method')
        for key in ('target_references', 'neanderthal_references', 'denisovan_references'):
            if len(set(g[key])) != len(g[key]):
                raise ValueError(f'Duplicate reference in gmm.{key}')
    for pair in c['plots']['contour_pairs']:
        if len(pair) != 2:
            raise ValueError('Contour requires exactly two references')
        used += pair
    if not set(used) <= set(ids):
        raise ValueError(f'Unknown references: {set(used) - set(ids)}')
    groups = {r['id']: r['group'] for r in refs}
    if any(groups[r] != 'neanderthal' for r in cl['neanderthal_references']):
        raise ValueError('Classification Neanderthal references must have group neanderthal')
    if any(groups[r] != 'denisovan' for r in cl['denisovan_references']):
        raise ValueError('Classification Denisovan references must have group denisovan')
    if g['enabled']:
        for key, group in [('target_references', 'denisovan'), ('denisovan_references', 'denisovan'), ('neanderthal_references', 'neanderthal')]:
            if any(groups[r] != group for r in g[key]):
                raise ValueError(f'Wrong group in gmm.{key}')
    c['project_root'] = str(ROOT)
    return c


def main():
    """Parse CLI arguments, save the effective config and invoke Snakemake."""
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['run'])
    p.add_argument('--config', default=str(ROOT / 'config/config.yaml'))
    for key in ('vcf', 'samples', 'outdir'):
        p.add_argument('--' + key)
    p.add_argument('--cores', type=int, default=4)
    p.add_argument('--memory-mb', type=int, default=32000)
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--target', choices=['all', 'matches', 'gmm', 'individual'], default='all')
    a = p.parse_args()
    if a.cores < 1 or a.memory_mb < 1:
        p.error('cores and memory-mb must be positive')
    c = resolve_config(a.config, vars(a))
    out = Path(c['outdir'])
    out.mkdir(parents=True, exist_ok=True)
    # Saving the fully resolved configuration makes each run reproducible and
    # lets Snakemake consume JSON without resolving paths a second time.
    cfg = out / 'resolved_config.json'
    content = json.dumps(c, indent=2) + '\n'
    if not cfg.exists() or cfg.read_text() != content:
        cfg.write_text(content)
    # Snakemake owns all parallelism. Individual scripts do not start their own
    # GNU Parallel pools, preventing accidental CPU oversubscription.
    cmd = ['snakemake', '--snakefile', str(ROOT / 'workflow/Snakefile'),
           '--directory', str(out), '--configfile', str(cfg), '--cores', str(a.cores),
           '--resources', f'mem_mb={a.memory_mb}', '--rerun-incomplete',
           # All rules use one dispatcher script. Its whole-file Snakemake
           # code hash would invalidate unrelated branches, so code changes
           # are represented by the stage-local ``params.code`` signatures.
           '--rerun-triggers', 'input', 'mtime', 'params', 'software-env',
           '--printshellcmds']
    if a.dry_run:
        cmd.append('--dry-run')
    cmd.append(a.target)
    return subprocess.call(cmd)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, FileNotFoundError) as e:
        sys.exit(str(e))
