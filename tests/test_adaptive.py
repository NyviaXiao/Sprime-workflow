import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('adaptive', ROOT/'workflow/scripts/run_adaptive.py')
adaptive = importlib.util.module_from_spec(spec); spec.loader.exec_module(adaptive)

class AdaptiveTests(unittest.TestCase):
    def test_af_and_core(self):
        self.assertAlmostEqual(adaptive.introgressed_af(.7, 1), .7)
        self.assertAlmostEqual(adaptive.introgressed_af(.7, 0), .3)
        rows = [{'introgressed_AF': x} for x in (.25, .32, .50, .65)]
        core, maximum = adaptive.core_variants(rows)
        self.assertEqual([r['introgressed_AF'] for r in core], [.50, .65])
        self.assertEqual(maximum, .65)

    def test_reference_thresholds_and_top2(self):
        self.assertFalse(adaptive.adaptive_pass(5, .90, 'denisovan'))
        self.assertFalse(adaptive.adaptive_pass(12, .40, 'neanderthal'))
        self.assertTrue(adaptive.adaptive_pass(11, .60, 'neanderthal'))
        rows = [
            {'population':'P','segment_id':'A','chromosome':'1','candidate_start':0,'core_mean_AF':.80,'pass_flag':1},
            {'population':'P','segment_id':'B','chromosome':'1','candidate_start':10,'core_mean_AF':.75,'pass_flag':1},
            {'population':'P','segment_id':'C','chromosome':'1','candidate_start':20,'core_mean_AF':.70,'pass_flag':1},
        ]
        selected = adaptive.rank_top2(rows)
        self.assertEqual([r['segment_id'] for r in selected], ['A','B'])
        self.assertEqual(rows[2]['selected_top2'], '0')

    def test_shared_overlap(self):
        rows = [
            {'chromosome':'1','start':0,'end':100,'population':'POP1'},
            {'chromosome':'1','start':20,'end':80,'population':'POP2'},
            {'chromosome':'1','start':90,'end':120,'population':'POP3'},
        ]
        got = adaptive.shared_intervals(rows)
        self.assertEqual([(x['start'],x['end'],x['populations']) for x in got], [(20,80,'POP1,POP2'),(90,100,'POP1,POP3')])

if __name__ == '__main__': unittest.main()
