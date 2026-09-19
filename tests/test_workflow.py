import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import json

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


class WorkflowTests(unittest.TestCase):
    def test_contour_pair_modes(self):
        source = (ROOT/'workflow/Snakefile').read_text().split('wildcard_constraints:')[0]
        refs = [{'id': x, 'group': 'neanderthal' if x.startswith('n') else 'denisovan'} for x in ['n1','n2','d1','d2','d3']]
        def pairs(mode, explicit):
            cfg = {'populations':['P'], 'chromosomes':['1'], 'archaic_references':refs,
                   'project_root':str(ROOT), 'git_commit':'test',
                   'plots':{'contour':{'enabled':True,'mode':mode,'pairs':explicit}}}
            namespace = {'config':cfg}
            exec(source, namespace)
            return namespace['PAIR_IDS']
        self.assertEqual(len(pairs('all_pairs', [])), 10)
        self.assertEqual(pairs('explicit_pairs', [['n1','d1'],['d1','d2']]), [('n1','d1'),('d1','d2')])

    def test_full_dag_with_placeholder_inputs(self):
        if importlib.util.find_spec('snakemake') is None:
            self.skipTest('Snakemake is not installed')
        from snakemake.api import SnakemakeApi
        from snakemake.settings.types import ResourceSettings, ConfigSettings, DAGSettings, StorageSettings
        from run import resolve_config
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            c = resolve_config(ROOT/'config/config.yaml', {})
            c['chromosomes'] = ['1', '2']
            c['outdir'] = str(tmp)

            for k in ('genetic_map', 'vcf'):
                c[k] = str(tmp/k)
                Path(c[k]).touch()

            for ref in c['archaic_references']:
                for k in ('vcf','mask'):
                    ref[k] = str(tmp/(ref['id']+k))
                    Path(ref[k]).touch()

            (tmp / 'resolved_config.json').write_text(
                json.dumps(c, indent=2) + '\n'
            )
                    
            with patch.dict(os.environ, {'XDG_CACHE_HOME':str(tmp/'cache')}):
                with SnakemakeApi() as api:
                    wf = api.workflow(ResourceSettings(cores=2), config_settings=ConfigSettings(config=c),
                                      storage_settings=StorageSettings(shared_fs_usage=set()),
                                      snakefile=ROOT/'workflow/Snakefile', workdir=tmp)
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output):
                        wf.dag(DAGSettings(targets={'all'})).printdag()
                    dag = output.getvalue()
                    for name in ('sprime','archaic_match','fit_gmm','individual_chr',
                                 'classification_summary','individual_summary','contour',
                                 'affinity_table','introgression_landscape','affinity_landscape','adaptive_chr','collect_adaptive',
                                 'provenance_base','report'):
                        self.assertIn(name,dag)

    def test_report_contour_order_is_canonical(self):
        import yaml
        from run import resolve_config
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            config = yaml.safe_load((ROOT/'config/config.yaml').read_text())
            config['samples'] = str((ROOT/'config/samples.tsv').resolve())
            config['report']['contour_pair'] = ['altai_nean', 'denisovan3']
            good = tmp/'good.yaml'; good.write_text(yaml.safe_dump(config))
            self.assertEqual(resolve_config(good, {})['report']['contour_pair'], ['altai_nean', 'denisovan3'])
            config['report']['contour_pair'] = ['denisovan3', 'altai_nean']
            bad = tmp/'bad.yaml'; bad.write_text(yaml.safe_dump(config))
            with self.assertRaises(ValueError):
                resolve_config(bad, {})

    def test_provenance_rule_declares_freshness_inputs(self):
        rules = (ROOT/'workflow/rules/new_modules.smk').read_text()
        snakefile = (ROOT/'workflow/Snakefile').read_text()
        self.assertIn("validation='qc/validation.tsv'", rules)
        self.assertIn("run_info='run_info.json'", rules)
        self.assertIn('git_commit=CURRENT_GIT_COMMIT', rules)
        self.assertIn("CURRENT_GIT_COMMIT = C['git_commit']", snakefile)

    def test_gmm_input_log_contains_every_output_wildcard(self):
        rules = (ROOT/'workflow/rules/gmm.smk').read_text()
        self.assertIn("log: 'logs/{population}/gmm_input.{ref}.log'", rules)
        self.assertNotIn("log: 'logs/{population}/gmm_input.log'", rules)

    def test_final_presentation_artifacts_and_contour_manifest_are_declared(self):
        snakefile = (ROOT/'workflow/Snakefile').read_text()
        provenance = (ROOT/'workflow/scripts/build_provenance.py').read_text()
        self.assertIn("'report/report.pdf'", snakefile)
        self.assertIn("'gmm/plots/{population}.{ref}.png'", snakefile)
        self.assertIn("enabled['contour']", provenance)

    def test_optional_modules_can_be_disabled_in_dag(self):
        if importlib.util.find_spec('snakemake') is None:
            self.skipTest('Snakemake is not installed')
        from snakemake.api import SnakemakeApi
        from snakemake.settings.types import ResourceSettings, ConfigSettings, DAGSettings, StorageSettings
        from run import resolve_config
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp); c = resolve_config(ROOT/'config/config.yaml', {})
            c['chromosomes'] = ['1']; c['outdir'] = str(tmp)
            for key in ('genetic_map', 'vcf'):
                c[key] = str(tmp/key); Path(c[key]).touch()
            for ref in c['archaic_references']:
                for key in ('vcf', 'mask'):
                    ref[key] = str(tmp/(ref['id']+key)); Path(ref[key]).touch()
            c['individual']['enabled'] = c['affinity']['enabled'] = c['adaptive']['enabled'] = c['report']['enabled'] = False
            c['plots']['contour']['enabled'] = False

            (tmp / 'resolved_config.json').write_text(
                json.dumps(c, indent=2) + '\n'
            )
            with patch.dict(os.environ, {'XDG_CACHE_HOME':str(tmp/'cache')}):
                with SnakemakeApi() as api:
                    wf = api.workflow(ResourceSettings(cores=2), config_settings=ConfigSettings(config=c), storage_settings=StorageSettings(shared_fs_usage=set()), snakefile=ROOT/'workflow/Snakefile', workdir=tmp)
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output): wf.dag(DAGSettings(targets={'all'})).printdag()
            dag = output.getvalue()
            for name in ('individual_chr', 'affinity_table', 'adaptive_chr', 'contour', 'report'):
                self.assertNotIn(name, dag)
            self.assertIn('provenance_base', dag)


if __name__ == '__main__':
    unittest.main()
