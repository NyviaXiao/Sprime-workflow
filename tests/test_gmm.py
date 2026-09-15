import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
AVAILABLE = all(importlib.util.find_spec(m) for m in ('numpy','scipy','sklearn','matplotlib'))


@unittest.skipUnless(AVAILABLE, 'Scientific dependencies are not installed')
class GMMTests(unittest.TestCase):
    def test_bimodal_and_insufficient(self):
        import numpy as np
        from core import write, read
        from run_gmm import analyze
        rng = np.random.default_rng(42)
        x = np.r_[rng.normal(.45,.025,80), rng.normal(.85,.025,80)]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            source = p/'input.tsv'
            s = SimpleNamespace(config={'populations':['POP'], 'gmm':dict(
                target_references=['deni'], method='legacy_chi_square', min_segments=10,
                n_init=5, seed=0, alpha=.05, bootstrap_replicates=5)},
                wildcards=SimpleNamespace(population='POP',ref='deni'), input=[str(source)],
                output=SimpleNamespace(plot=str(p/'plot.png'), selection=str(p/'selection.tsv'), components=str(p/'components.tsv')))
            write(source,['deni_match_rate'],[{'deni_match_rate':v} for v in x])
            analyze(s)
            result = read(s.output.selection)[0]
            self.assertEqual(result['status'],'ok')
            self.assertEqual(result['selected_components'],'2')
            self.assertEqual(result['second_comparison'],'2_vs_3')
            for field in ('aic_1','aic_2','aic_3','bic_1','bic_2','bic_3'):
                self.assertNotEqual(result[field], 'NA')
                self.assertTrue(float(result[field]) == float(result[field]))
            self.assertEqual(len(read(s.output.components)),2)
            self.assertGreater(Path(s.output.plot).stat().st_size,1000)
            write(source,['deni_match_rate'],[])
            analyze(s)
            self.assertEqual(read(s.output.selection)[0]['status'],'insufficient_data')
            self.assertEqual(read(s.output.components),[])

    def test_bootstrap_probability(self):
        import numpy as np
        from run_gmm import fit, compare
        x = np.random.default_rng(2).normal(.5,.1,40)
        cfg = dict(n_init=2, method='parametric_bootstrap', bootstrap_replicates=5)
        models = {k:fit(x,k,cfg,0) for k in (1,2)}
        lr,p = compare(x,1,2,models,cfg,3)
        self.assertGreaterEqual(lr,0)
        self.assertTrue(1/6 <= p <= 1)


if __name__ == '__main__':
    unittest.main()
