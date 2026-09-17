import importlib.util
import json
import pathlib
import tempfile
import unittest
from unittest.mock import patch

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

if __name__ == '__main__': unittest.main()
