#!/usr/bin/env python3
"""Write small, deterministic run and artifact manifests."""
import csv, hashlib, json, os, platform, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

def _dump(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')

def _git_commit(root):
    try:
        return subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return 'unavailable'

def build_base(config_path, outputs, root):
    cfg = json.loads(Path(config_path).read_text())
    refs = [r['id'] for r in cfg['archaic_references']]
    enabled = {k: bool(cfg.get(k, {}).get('enabled', False)) for k in ('gmm','individual','landscape','affinity','adaptive','report')}
    _dump(outputs['resolved'], cfg)
    _dump(outputs['run'], {'run_time': datetime.now(timezone.utc).isoformat(), 'genome_build': cfg['genome_build'], 'populations': cfg['populations'], 'chromosomes': cfg['chromosomes'], 'archaic_references': refs, 'git_commit': _git_commit(root), 'enabled_modules': enabled})
    versions = [('python', platform.python_version()), ('snakemake', 'workflow-managed'), ('bcftools', 'workflow-managed'), ('java', 'workflow-managed'), ('R', 'workflow-managed')]
    try:
        import numpy, scipy, sklearn, pysam
        versions += [('numpy', numpy.__version__), ('scipy', scipy.__version__), ('scikit-learn', sklearn.__version__), ('pysam', pysam.__version__)]
    except ImportError:
        pass
    Path(outputs['software']).parent.mkdir(parents=True, exist_ok=True)
    with open(outputs['software'],'w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['software','version']); w.writerows(versions)
    paths = [('vcf', cfg['vcf']), ('samples', cfg['samples']), ('genetic_map', cfg['genetic_map']), ('sprime_jar', cfg['sprime_jar']), ('map_arch', cfg['map_arch'])]
    paths += [(f"{r['id']}_vcf", r['vcf']) for r in cfg['archaic_references']] + [(f"{r['id']}_mask", r['mask']) for r in cfg['archaic_references']]
    with open(outputs['inputs'],'w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['input_name','path','size','mtime'])
        for name, templ in paths:
            path = str(templ).replace('{chrom}', cfg['chromosomes'][0]); p=Path(path)
            if p.exists(): w.writerow([name,path,p.stat().st_size,p.stat().st_mtime])

def build_outputs(paths, output):
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['path','size','sha256'])
        for path in sorted(map(str, paths)):
            p=Path(path)
            if not p.is_file(): raise FileNotFoundError(path)
            h=hashlib.sha256();
            with open(p,'rb') as data:
                for block in iter(lambda:data.read(1024*1024),b''): h.update(block)
            w.writerow([path,p.stat().st_size,h.hexdigest()])

if 'snakemake' in globals():
    if snakemake.params.mode == 'base':
        build_base(snakemake.input.config, {'run':snakemake.output.run,'resolved':snakemake.output.resolved,'software':snakemake.output.software,'inputs':snakemake.output.inputs}, snakemake.params.root)
    else:
        build_outputs(snakemake.input.artifacts, snakemake.output[0])
