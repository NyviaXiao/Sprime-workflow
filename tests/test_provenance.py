import importlib.util
import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch
import types
import sys

import run

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('provenance', ROOT/'workflow/scripts/build_provenance.py')
provenance = importlib.util.module_from_spec(spec); spec.loader.exec_module(provenance)

class ProvenanceTests(unittest.TestCase):
    def test_declared_inputs_expand_every_chromosome(self):
        cfg={'samples':'samples.tsv','genetic_map':'map','sprime_jar':'jar','map_arch':'arch','vcf':'modern.{chrom}.vcf.gz','chromosomes':['1','2'],'archaic_references':[{'id':'n1','vcf':'n1.{chrom}.vcf.gz','mask':'n1.{chrom}.bed.gz'}]}
        rows=dict(provenance._declared_inputs(cfg))
        self.assertEqual(rows['vcf_chr1'],'modern.1.vcf.gz'); self.assertEqual(rows['vcf_chr2'],'modern.2.vcf.gz')
        self.assertEqual(rows['n1_mask_chr2'],'n1.2.bed.gz')

    def test_run_info_summary_and_git_dirty_state(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory) / 'run_info.json'
            path.write_text(json.dumps({'time_utc':'now','versions':{'R':'R version'},'sha256':{'samples':'abc'}}))
            summary = provenance._run_info_summary(path)
        self.assertEqual(summary['tool_versions']['R'], 'R version')
        clean = ["abc\n", "", b'', b'']
        with patch('run.subprocess.check_output', side_effect=clean):
            metadata = run.git_metadata(pathlib.Path('/repo'))
        self.assertFalse(metadata['git_dirty'])
        dirty = ["abc\n", " M workflow.py\n", b'one', b'two']
        with patch('run.subprocess.check_output', side_effect=dirty):
            metadata = run.git_metadata(pathlib.Path('/repo'))
        self.assertTrue(metadata['git_dirty'])
        self.assertEqual(metadata['git_diff_sha256'], run.hashlib.sha256(b'onetwo').hexdigest())

    def test_build_base_writes_all_provenance_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            paths = {name: root / name for name in ('samples.tsv', 'map.txt', 'sprime.jar', 'map_arch', 'modern.vcf', 'n1.vcf', 'n1.bed')}
            for path in paths.values(): path.write_text('x')
            cfg = {'samples':str(paths['samples.tsv']), 'genetic_map':str(paths['map.txt']),
                   'sprime_jar':str(paths['sprime.jar']), 'map_arch':str(paths['map_arch']),
                   'vcf':str(paths['modern.vcf']), 'chromosomes':['1'],
                   'archaic_references':[{'id':'n1','vcf':str(paths['n1.vcf']),'mask':str(paths['n1.bed'])}],
                   'genome_build':'GRCh37', 'populations':['P'], 'plots':{'contour':{'enabled':True}},
                   'git_dirty':False, 'git_diff_sha256':'digest'}
            config, run_info = root/'resolved_config.json', root/'run_info.json'
            config.write_text(json.dumps(cfg)); run_info.write_text(json.dumps({'time_utc':'now','versions':{},'sha256':{}}))
            outputs = {'run':root/'run_manifest.json', 'resolved':root/'provenance_config.json',
                       'software':root/'software.tsv', 'inputs':root/'inputs.tsv'}
            fake_modules = {name:types.SimpleNamespace(__version__='test') for name in ('numpy','pysam','scipy','sklearn')}
            with patch.object(provenance, '_command_version', return_value='version'), patch.dict(sys.modules, fake_modules):
                provenance.build_base(config, run_info, outputs, 'commit')
            for path in outputs.values():
                self.assertTrue(path.exists())
                self.assertTrue(path.read_text())
            self.assertEqual(json.loads(outputs['run'].read_text())['git_commit'], 'commit')

if __name__ == '__main__': unittest.main()
