import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


@unittest.skipUnless(importlib.util.find_spec('snakemake'), 'Snakemake is not installed')
class WorkflowTests(unittest.TestCase):
    def test_full_dag_with_placeholder_inputs(self):
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
            with patch.dict(os.environ, {'XDG_CACHE_HOME':str(tmp/'cache')}):
                with SnakemakeApi() as api:
                    wf = api.workflow(ResourceSettings(cores=2), config_settings=ConfigSettings(config=c),
                                      storage_settings=StorageSettings(shared_fs_usage=set()),
                                      snakefile=ROOT/'workflow/Snakefile', workdir=tmp)
                    output = io.StringIO()
                    with contextlib.redirect_stdout(output):
                        wf.dag(DAGSettings(targets={'all'})).printdag()
                    dag = output.getvalue()
                    for name in ('sprime','archaic_match','fit_gmm','individual_chr','contour'):
                        self.assertIn(name,dag)


if __name__ == '__main__':
    unittest.main()
