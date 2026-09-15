"""Execute one Snakemake job and send all messages to that job's log.

Every rule passes a short ``params.stage`` value to this dispatcher. Keeping
external-command handling here gives bcftools, Java, map_arch and R identical
error handling and logging, while scientific helpers remain in focused files.
Commands are argv lists, so paths and sample names are never shell-expanded.
"""
import contextlib
import os
import subprocess
import sys
import traceback
from collections import defaultdict
from pathlib import Path

# Snakemake executes a generated script; derive module location from config.
sys.path.insert(0, str(Path(snakemake.config['project_root']) / 'workflow/scripts'))
from core import (BASE, CALLS, classify, classification_summary_rows, collect,
                  gmm_target_pass, individual_summary_rows, read, score_rows,
                  summarize, summary_fields, write)
from validate_inputs import validate
from call_individual_tracts import individuals


def command(args, output=None):
    """Run an external command and raise immediately on a non-zero exit code.

    When ``output`` is given, only stdout is redirected there. Diagnostic
    stderr remains in the Snakemake log surrounding this dispatcher.
    """
    print('COMMAND', repr([str(a) for a in args]), flush=True)
    if output:
        with open(output, 'w') as f:
            subprocess.run(list(map(str, args)), check=True, stdout=f, stderr=sys.stderr)
    else:
            subprocess.run(list(map(str, args)), check=True, stdout=sys.stdout, stderr=sys.stderr)


def as_paths(value):
    """Normalize a named Snakemake input to a list without changing paths."""
    return [value] if isinstance(value, (str, Path)) else list(value)


def main(s):
    """Dispatch the current Snakemake job to its implementation stage."""
    c, stage = s.config, s.params.stage
    refs = {r['id']: r for r in c['archaic_references']}

    # Tool preparation and input validation.
    if stage == 'compile':
        command(['make', '-B', '-C', Path(s.input.source).parent])
    elif stage == 'validate':
        validate(s, c)

    # Prepare one VCF containing target-population and shared outgroup samples.
    elif stage == 'samples':
        rows = read(s.input.samples)
        outgroup = [r['sample_id'] for r in rows if r['role'] == 'outgroup']
        target = [r['sample_id'] for r in rows
                  if r['role'] == 'target'
                  and r['population'] == s.wildcards.population]
        Path(s.output.samples).write_text('\n'.join(target + outgroup) + '\n')
        Path(s.output.outgroup).write_text('\n'.join(outgroup) + '\n')
    elif stage == 'subset':
        import tempfile
        # Filtering through temporary BCF avoids repeatedly compressing data.
        # Final files contain biallelic SNPs, GT only, and at least one ALT.
        with tempfile.TemporaryDirectory(dir=Path(s.output.vcf).parent) as tmp:
            subset, filtered = Path(tmp)/'subset.bcf', Path(tmp)/'filtered.bcf'
            command(['bcftools', 'view', '-t', s.wildcards.chrom,
                     '-S', s.input.samples, '-Ob', '-o', subset, s.input.vcf])
            command(['bcftools', 'view', '-c', '1', '-m2', '-M2', '-v',
                     'snps', '-Ob', '-o', filtered, subset])
            command(['bcftools', 'annotate', '-x', 'INFO,^FORMAT/GT', '-Oz',
                     '-o', s.output.vcf, filtered])
        command(['bcftools', 'index', '--csi', s.output.vcf])
    elif stage == 'concat':
        # Inputs arrive in the chromosome order declared in config.yaml.
        command(['bcftools', 'concat', '-Oz', '-o', s.output.vcf, *s.input.vcfs])

    # Detect SPrime candidate segments, then annotate each ancient reference
    # independently. Independent annotation permits parallel/retry per reference.
    elif stage == 'sprime':
        prefix = str(s.output.score).removesuffix('.score')
        # Leave 20% of the rule's memory allocation for Java/native overhead.
        heap = max(256, int(s.resources.mem_mb * .8))
        command(['java', f'-Xmx{heap}m', '-jar', s.input.jar, f'gt={s.input.vcf}',
                 f'outgroup={s.input.outgroup}', f'map={s.input.map}', f'out={prefix}',
                 f'chrom={s.wildcards.chrom}', f'minscore={c["sprime_minscore"]}'])
        if not Path(s.output.score).exists():
            raise RuntimeError('SPrime did not produce a score file; inspect log')
    elif stage == 'match':
        # map_arch cannot process a header-only score file, so create the valid
        # empty annotated output directly when SPrime found no markers.
        if next(score_rows(s.input.score), None) is None:
            with open(s.input.score) as f:
                header = f.readline().strip().split()
            Path(s.output.score).write_text('\t'.join(header + [s.params.tag]) + '\n')
        else:
            command([s.input.binary, '--kp', '--sep', '\\t', '--tag', s.params.tag,
                     '--mskbed', s.input.mask, '--vcf', s.input.vcf, '--score', s.input.score], s.output.score)
    # Build normal summary rates and a separate pre_data-compatible GMM table.
    # The latter applies SCORE > threshold before recounting callable alleles.
    elif stage in ('summarize', 'gmm_base'):
        fields = summary_fields(refs)
        minimum = c['summary_min_callable'] if stage == 'summarize' else 1
        threshold = None if stage == 'summarize' else c['gmm']['score_gt']
        write(s.output[0], fields, summarize(s.input, s.wildcards.population, refs, minimum, threshold))
    elif stage == 'collect_tables':
        # Wide output is convenient for users; long output is tidy and easier
        # to aggregate or plot for an arbitrary number of ancient references.
        collect(s.input, s.output.wide)
        long_fields = BASE + ['archaic_id', 'archaic_group', 'matched',
                              'mismatch', 'callable', 'notcomp', 'match_rate']
        totals = defaultdict(lambda: [0, 0, 0])

        def emit():
            for path in s.input:
                for row in read(path):
                    for rid, ref in refs.items():
                        extra = {k: row[f'{rid}_{k}'] for k in
                                 ('matched', 'mismatch', 'callable',
                                  'notcomp', 'match_rate')}
                        bucket = totals[row['population'], rid]
                        bucket[0] += 1
                        bucket[1] += int(extra['matched'])
                        bucket[2] += int(extra['callable'])
                        yield {**row, **extra, 'archaic_id': rid, 'archaic_group': ref['group']}
        write(s.output.long, long_fields, emit())
        pooled_fields = ['population', 'archaic_id', 'n_segments', 'matched',
                         'callable', 'pooled_match_rate']
        pooled_rows = [dict(zip(pooled_fields,
                                [*key, *value,
                                 value[1]/value[2] if value[2] else 'NA']))
                       for key, value in sorted(totals.items())]
        write(s.output.summary, pooled_fields, pooled_rows)

    # Segment classification and GMM model selection.
    elif stage == 'classify':
        rows = read(s.input[0])
        write(s.output[0], summary_fields(refs) + ['archaic_class'],
              ({**r, 'archaic_class': classify(r, c['classification'])} for r in rows))
    elif stage == 'classification_summary':
        # Count the already classified segments and sum their displayed marker
        # span. Emit all three classes for every population, including zeros.
        fields = ['population', 'archaic_class', 'n_segments', 'genome_coverage_mb']
        class_rows = {population: read(path)
                      for population, path in zip(c['populations'], as_paths(s.input.classification))}
        length_rows = {population: read(path)
                       for population, path in zip(c['populations'], as_paths(s.input.tables))}
        write(s.output[0], fields,
              classification_summary_rows(class_rows, length_rows, c['populations']))
    elif stage == 'gmm_input':
        target = s.wildcards.ref
        passing_rows = (row for row in read(s.input[0])
                        if gmm_target_pass(row, c['gmm'], target))
        write(s.output[0], summary_fields(refs), passing_rows)
    elif stage == 'gmm':
        from run_gmm import analyze
        analyze(s)
    elif stage == 'collect_gmm':
        collect(s.input.selection, s.output.selection)
        collect(s.input.components, s.output.components)
    # Individual calls are computed per chromosome and then split into the
    # three mutually exclusive archaic classes at population level.
    elif stage == 'individual':
        individuals(s, c)
    elif stage == 'individual_pop':
        collect(s.input, s.output.all)
        rows = read(s.output.all)
        for output in list(s.output)[1:]:
            label = Path(output).name.split('.')[0]
            write(output, CALLS, (r for r in rows if r['archaic_class'] == label))
    elif stage == 'individual_summary':
        rows_by_population = {population: read(path)
                              for population, path in zip(c['populations'], as_paths(s.input.calls))}
        samples = read(s.input.samples)
        fields = ['population', 'archaic_class', 'per_individual_introgressed_mb',
                  'n_introgressed_individuals']
        write(s.output[0], fields,
              individual_summary_rows(rows_by_population, samples, c['populations']))
    elif stage == 'collect':
        collect(s.input, s.output[0])
    elif stage == 'contour':
        a, b = s.wildcards.pair.split('__')
        command(['Rscript', s.input.script, s.input.table, s.output[0],
                 a, b, s.wildcards.population])
    else:
        raise ValueError(stage)


# A Snakemake ``script:`` receives the global ``snakemake`` object. Prepare
# output/log directories here so each stage can focus on its own operation.
for output in snakemake.output:
    Path(output).parent.mkdir(parents=True, exist_ok=True)
Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.log[0], 'w') as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
    # Prevent NumPy/BLAS from creating hidden thread pools that exceed the
    # resource allocation controlled by Snakemake.
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[key] = '1'
    try:
        main(snakemake)
    except Exception:
        traceback.print_exc()
        raise
