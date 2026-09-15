import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


class SignatureTests(unittest.TestCase):
    def test_stage_signatures_are_not_global(self):
        from run import resolve_config

        namespace = {'config': resolve_config(ROOT / 'config/config.yaml', {})}
        source = (ROOT / 'workflow/Snakefile').read_text(encoding='utf-8')
        # Evaluate only the ordinary Python preamble; the remainder uses the
        # Snakemake DSL and is intentionally not imported in this unit test.
        exec(source.split('wildcard_constraints:')[0], namespace)
        signatures = namespace['SIG']
        self.assertNotEqual(signatures['gmm'], signatures['classify'])
        self.assertNotEqual(signatures['gmm'], signatures['individual'])
        self.assertNotEqual(signatures['classify'], signatures['individual'])
        self.assertEqual(len(signatures), 20)


if __name__ == '__main__':
    unittest.main()
