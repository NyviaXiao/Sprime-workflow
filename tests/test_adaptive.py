import csv
import importlib.util
import pathlib
import tempfile
import sys
import types
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
            self.assertEqual(top[0]['adaptive_group'], 'neanderthal')

    def test_group_specific_top2_and_shared_overlaps(self):
        rows = [
            {'population':'P','chromosome':'1','segment_id':'bad','candidate_start':'0','candidate_end':'100','core_mean_AF':'0.99','neanderthal_pass':'0','denisovan_pass':'0'},
            {'population':'P','chromosome':'1','segment_id':'n1','candidate_start':'200','candidate_end':'300','core_mean_AF':'0.80','neanderthal_pass':'1','denisovan_pass':'0'},
            {'population':'P','chromosome':'1','segment_id':'n2','candidate_start':'400','candidate_end':'500','core_mean_AF':'0.70','neanderthal_pass':'1','denisovan_pass':'0'},
            {'population':'P','chromosome':'1','segment_id':'n3','candidate_start':'600','candidate_end':'700','core_mean_AF':'0.60','neanderthal_pass':'1','denisovan_pass':'0'},
            {'population':'P','chromosome':'1','segment_id':'dual','candidate_start':'800','candidate_end':'900','core_mean_AF':'0.95','neanderthal_pass':'1','denisovan_pass':'1'},
            {'population':'P','chromosome':'1','segment_id':'d1','candidate_start':'1000','candidate_end':'1100','core_mean_AF':'0.85','neanderthal_pass':'0','denisovan_pass':'1'},
        ]
        self.assertEqual([row['segment_id'] for row in collector.rank_top2(rows, 'neanderthal')], ['dual', 'n1'])
        self.assertEqual([row['segment_id'] for row in collector.rank_top2(rows, 'denisovan')], ['dual', 'd1'])
        self.assertEqual(collector.rank_top2(rows, 'denisovan')[0]['adaptive_group'], 'denisovan')
        overlap = collector.shared_intervals([
            {'adaptive_group':'neanderthal','population':'P1','chromosome':'1','candidate_start':0,'candidate_end':100},
            {'adaptive_group':'neanderthal','population':'P2','chromosome':'1','candidate_start':20,'candidate_end':80},
            {'adaptive_group':'neanderthal','population':'P3','chromosome':'1','candidate_start':90,'candidate_end':120},
            {'adaptive_group':'denisovan','population':'P4','chromosome':'1','candidate_start':20,'candidate_end':80},
        ])
        self.assertEqual([(row['start'], row['end'], row['populations']) for row in overlap], [(20,80,'P1,P2'), (90,100,'P1,P3')])
        self.assertTrue(all(row['adaptive_group'] == 'neanderthal' for row in overlap))

    def test_collect_writes_two_ranked_groups_and_allows_dual_pass(self):
        fields = ['population','chromosome','segment_id','candidate_start','candidate_end',
                  'core_mean_AF','neanderthal_pass','denisovan_pass','pass_flag']
        rows = [
            {'population':'P','chromosome':'1','segment_id':'dual','candidate_start':0,'candidate_end':1,'core_mean_AF':.9,'neanderthal_pass':1,'denisovan_pass':1,'pass_flag':1},
            {'population':'P','chromosome':'1','segment_id':'n','candidate_start':2,'candidate_end':3,'core_mean_AF':.8,'neanderthal_pass':1,'denisovan_pass':0,'pass_flag':1},
            {'population':'P','chromosome':'1','segment_id':'d','candidate_start':4,'candidate_end':5,'core_mean_AF':.7,'neanderthal_pass':0,'denisovan_pass':1,'pass_flag':1},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            variants, core, segments = tmp/'adaptive.1.variants.tsv.gz', tmp/'adaptive.1.core.tsv.gz', tmp/'adaptive.1.segments.tsv'
            adaptive._write(variants, ['population'], []); adaptive._write(core, ['population'], [])
            adaptive._write(segments, fields, rows)
            outputs = {name:tmp/f'{name}.tsv' for name in ('variants','core','segments','top2','shared')}
            collector.collect([variants, core, segments], outputs)
            with open(outputs['top2'], newline='') as handle:
                top = list(csv.DictReader(handle, delimiter='\t'))
            with open(outputs['segments'], newline='') as handle:
                segment_fields = csv.DictReader(handle, delimiter='\t').fieldnames
            top_fields = list(top[0])
        self.assertEqual([(row['adaptive_group'], row['segment_id'], row['adaptive_rank']) for row in top],
                         [('neanderthal','dual','1'), ('neanderthal','n','2'),
                          ('denisovan','dual','1'), ('denisovan','d','2')])
        for name in ('adaptive_group', 'adaptive_rank', 'selected_top2'):
            self.assertNotIn(name, segment_fields)
            self.assertIn(name, top_fields)

    def test_alt_frequency_keeps_only_sprime_keys(self):
        class Record:
            def __init__(self, pos):
                self.pos, self.ref, self.alts = pos, 'A', ('G',)
                self.samples = {'sample': {'GT': (0, 1)}}
        class VariantFile:
            def __init__(self, path): pass
            def __enter__(self): return self
            def __exit__(self, *args): pass
            def fetch(self, chrom): return [Record(10), Record(20)]
        original = sys.modules.get('pysam')
        sys.modules['pysam'] = types.SimpleNamespace(VariantFile=VariantFile)
        try:
            got = adaptive._alt_frequencies('unused.vcf', '1', ['sample'], {(10, 'A', 'G')})
        finally:
            if original is None: del sys.modules['pysam']
            else: sys.modules['pysam'] = original
        self.assertEqual(set(got), {(10, 'A', 'G')})


if __name__ == '__main__': unittest.main()
