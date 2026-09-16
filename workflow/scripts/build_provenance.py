#!/usr/bin/env python3
"""Write reproducible run/input/output manifests from declared workflow files."""
import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def _dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def _command_version(command):
    completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, check=True)
    lines = (completed.stdout + '\n' + completed.stderr).splitlines()
    return next((line.strip() for line in lines if line.strip()), '')


def _git_commit(root):
    try:
        return subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unavailable'


def _declared_inputs(cfg):
    items = [('samples', cfg['samples']), ('genetic_map', cfg['genetic_map']),
             ('sprime_jar', cfg['sprime_jar']), ('map_arch', cfg['map_arch'])]
    templates = [('vcf', cfg['vcf'])]
    for ref in cfg['archaic_references']:
        templates.extend([(f"{ref['id']}_vcf", ref['vcf']),
                          (f"{ref['id']}_mask", ref['mask'])])
    for name, template in templates:
        if '{chrom}' in str(template):
            items.extend((f'{name}_chr{chrom}', str(template).format(chrom=chrom))
                         for chrom in cfg['chromosomes'])
        else:
            items.append((name, template))
    return items


def build_base(config_path, outputs, root):
    cfg = json.loads(Path(config_path).read_text())
    enabled = {name: bool(cfg.get(name, {}).get('enabled', False))
               for name in ('gmm', 'individual', 'landscape', 'affinity', 'adaptive', 'report')}
    _dump(outputs['resolved'], cfg)
    _dump(outputs['run'], {
        'run_time': datetime.now(timezone.utc).isoformat(),
        'genome_build': cfg['genome_build'], 'populations': cfg['populations'],
        'chromosomes': cfg['chromosomes'],
        'archaic_references': [ref['id'] for ref in cfg['archaic_references']],
        'git_commit': _git_commit(root), 'enabled_modules': enabled,
    })
    versions = [('python', _command_version([sys.executable, '--version'])),
                ('snakemake', _command_version(['snakemake', '--version'])),
                ('bcftools', _command_version(['bcftools', '--version'])),
                ('java', _command_version(['java', '-version'])),
                ('R', _command_version(['R', '--version']))]
    import numpy, pysam, scipy, sklearn
    versions.extend([('numpy', numpy.__version__), ('scipy', scipy.__version__),
                     ('scikit-learn', sklearn.__version__), ('pysam', pysam.__version__)])
    with open(outputs['software'], 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['software', 'version']); writer.writerows(versions)
    with open(outputs['inputs'], 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['input_name', 'path', 'size', 'mtime'])
        for name, path in _declared_inputs(cfg):
            file_path = Path(path)
            if not file_path.is_file():
                raise FileNotFoundError(f'Declared input missing: {file_path}')
            stat = file_path.stat()
            writer.writerow([name, file_path, stat.st_size, stat.st_mtime])


def build_outputs(paths, output):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['path', 'size', 'sha256'])
        for path in sorted(map(str, paths)):
            file_path = Path(path)
            if not file_path.is_file():
                raise FileNotFoundError(path)
            digest = hashlib.sha256()
            with open(file_path, 'rb') as data:
                for block in iter(lambda: data.read(1024 * 1024), b''):
                    digest.update(block)
            writer.writerow([file_path, file_path.stat().st_size, digest.hexdigest()])


if 'snakemake' in globals():
    if snakemake.params.mode == 'base':
        build_base(snakemake.input.config, {
            'run': snakemake.output.run, 'resolved': snakemake.output.resolved,
            'software': snakemake.output.software, 'inputs': snakemake.output.inputs,
        }, snakemake.params.root)
    else:
        build_outputs(snakemake.input.artifacts, snakemake.output[0])
