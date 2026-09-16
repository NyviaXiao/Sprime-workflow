import csv
import importlib.util
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).parents[1]
def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/f'workflow/scripts/{name}.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module
adaptive = load('run_adaptive')
collector = load('collect_adaptive')


class AdaptiveTests(unittest.TestCase):
    def test_af_and_core(self):
        self.assertAlmostEqual(adaptive.introgressed_af(.7, 1), .7)
        self.assertAlmostEqual(adaptive.introgressed_af(.7, 0), .3)
        core, maximum = adaptive.core_variants([{'introgressed_AF': x} for x in (.25, .32, .50, .65)], .30, .20)
        self.assertEqual([row['introgressed_AF'] for row in core], [.50, .65])
        self.assertEqual(maximum, .65)

    def test_chromosome_to_collect_integration(self):
        cfg = {'min_allele_frequency': .30, 'max_frequency_drop': .20,
               'min_callable_sites': 10, 'neanderthal_match_threshold': .50,
               'denisovan_match_threshold': .40}
        refs = ['n1', 'n2', 'd1', 'd2']
        scores = [{'POS': str(101+i), 'REF': 'A', 'ALT': 'G', 'ALLELE': '1', 'SEGMENT': 's1'} for i in range(10)]
        af = {(101+i, 'A', 'G'): .80 for i in range(10)}
        states = {ref: {(101+i, 'A', 'G'): ('match' if ref == 'n2' else 'mismatch') for i in range(10)} for ref in refs}
        markers, cores, segments = adaptive.build_chromosome_rows('P', '1', scores, af, states, cfg, refs, {'n1':'neanderthal','n2':'neanderthal','d1':'denisovan','d2':'denisovan'}, ['n1','n2'], ['d1','d2'])
        self.assertEqual(len(markers), 10); self.assertEqual(len(cores), 10)
        self.assertEqual(segments[0]['candidate_start'], 100)
        self.assertEqual(segments[0]['candidate_end'], 110)
        self.assertEqual(segments[0]['pass_flag'], 1)
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            paths = [tmp/'adaptive.1.variants.tsv.gz', tmp/'adaptive.1.core.tsv.gz', tmp/'adaptive.1.segments.tsv']
            adaptive._write(paths[0], list(markers[0]), markers); adaptive._write(paths[1], list(cores[0]), cores); adaptive._write(paths[2], list(segments[0]), segments)
            outputs = {'variants':tmp/'variants.tsv.gz','core':tmp/'core.tsv.gz','segments':tmp/'segments.tsv','top2':tmp/'top2.tsv','shared':tmp/'shared.tsv'}
            collector.collect(paths, outputs)
            for path in outputs.values(): self.assertTrue(path.exists())
            with open(outputs['top2'], newline='') as handle: top = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(top[0]['segment_id'], 's1')

    def test_string_zero_never_ranks_or_overlaps(self):
        rows = [{'population':'P','chromosome':'1','segment_id':'bad','candidate_start':'0','candidate_end':'100','core_mean_AF':'0.99','pass_flag':'0'}, {'population':'P','chromosome':'1','segment_id':'good','candidate_start':'200','candidate_end':'300','core_mean_AF':'0.80','pass_flag':'1'}]
        self.assertEqual([row['segment_id'] for row in collector.rank_top2(rows)], ['good'])
        overlap = collector.shared_intervals([{'population':'P1','chromosome':'1','candidate_start':0,'candidate_end':100}, {'population':'P2','chromosome':'1','candidate_start':20,'candidate_end':80}, {'population':'P3','chromosome':'1','candidate_start':90,'candidate_end':120}])
        self.assertEqual([(row['start'], row['end'], row['populations']) for row in overlap], [(20,80,'P1,P2'), (90,100,'P1,P3')])


if __name__ == '__main__': unittest.main()
