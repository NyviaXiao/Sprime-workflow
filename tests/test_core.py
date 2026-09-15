import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from core import (classify, classification_summary_rows, collect, gmm_pass,
                  gmm_target_pass, individual_summary_rows, read, runs,
                  summarize, summary_fields, write)


class CoreTests(unittest.TestCase):
    def test_runs_boundaries_missing_and_singletons(self):
        self.assertEqual(list(runs([(1,True),(9,True),(10,None),(20,True),(30,True)])), [(1,9,2),(20,30,2)])
        self.assertEqual(list(runs([(1,True),(2,False),(3,True)])), [])
        self.assertEqual(list(runs([(10,True)], 1)), [(10,10,1)])

    def test_summary_callable_score_and_empty(self):
        refs = {'nean': {'tag': 'Nean'}, 'deni': {'tag': 'Deni'}}
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for rid, tag in [('nean','Nean'), ('deni','Deni')]:
                path = Path(tmp)/rid/'chr1.mscore'
                write(path, ['CHROM','POS','ID','REF','ALT','SEGMENT','ALLELE','SCORE',tag], [
                    dict(zip(['CHROM','POS','ID','REF','ALT','SEGMENT','ALLELE','SCORE',tag],
                             ['1',pos,'.','A','G','7',1,score,state]))
                    for pos, score, state in [(10,150000,'match'),(20,160000,'mismatch'),(30,160000,'notcomp')]])
                paths.append(path)
            row = summarize(paths, 'POP', refs, 2)[0]
            self.assertEqual(row['nean_callable'], 2)
            self.assertEqual(row['nean_match_rate'], .5)
            self.assertEqual(row['marker_count'], 3)
            filtered = summarize(paths, 'POP', refs, 1, 150000)[0]
            self.assertEqual(filtered['nean_match_rate'], 0)
            self.assertEqual(filtered['marker_count'], 2)
            self.assertEqual(summarize(paths, 'POP', refs, 1, 200000), [])

    def test_gmm_no_length_filter_and_strict_boundaries(self):
        cfg = dict(neanderthal_references=['n'], denisovan_references=['d','d2'], target_references=['d'],
                   min_callable=30, neanderthal_rate_lt=.3, denisovan_rate_gt=.3)
        row = dict(n_callable=30, d_callable=30, d2_callable=30, n_match_rate=.2,
                   d_match_rate=.5, d2_match_rate=.1, length_bp=1)
        self.assertTrue(gmm_pass(row,cfg))
        self.assertFalse(gmm_pass({**row,'n_match_rate':.3},cfg))
        self.assertFalse(gmm_pass({**row,'d_callable':29},cfg))
        self.assertFalse(gmm_pass({**row,'d_callable':29,'d2_callable':29,'d2_match_rate':.9},cfg))
        # A non-callable Deni reference is excluded; the callable Deni remains
        # sufficient for candidate eligibility.
        self.assertTrue(gmm_pass({**row,'d2_match_rate':'NA'},cfg))
        self.assertTrue(gmm_pass({**row,'d2_callable':10,'d2_match_rate':.9},cfg))
        self.assertFalse(gmm_pass({**row,'d_match_rate':'NA','d2_match_rate':'NA'},cfg))
        self.assertFalse(gmm_pass({**row,'n_callable':29},cfg))
        self.assertTrue(gmm_target_pass({**row,'d2_callable':10},cfg,'d'))
        self.assertFalse(gmm_target_pass({**row,'d2_callable':10},cfg,'d2'))

    def test_multi_reference_classification_and_boundaries(self):
        cfg = dict(neanderthal_references=['n1','n2'], denisovan_references=['d1','d2'], neanderthal_gt=.6,
                   denisovan_lt_for_neanderthal=.4, denisovan_gt=.3, neanderthal_lt_for_denisovan=.3)
        self.assertEqual(classify({'n1_match_rate':.2,'n2_match_rate':.75,
                                   'd1_match_rate':.1,'d2_match_rate':.2},cfg), 'neanderthal')
        self.assertEqual(classify({'n1_match_rate':.1,'n2_match_rate':.2,
                                   'd1_match_rate':.2,'d2_match_rate':.8},cfg), 'denisovan')
        self.assertEqual(classify({'n1_match_rate':.8,'n2_match_rate':.2,
                                   'd1_match_rate':.7,'d2_match_rate':.1},cfg), 'ambiguous')
        for n,d in [(.6,.1),(.8,.4),(.1,.3),(.3,.8)]:
            self.assertEqual(classify({'n1_match_rate':n,'n2_match_rate':'NA',
                                       'd1_match_rate':d,'d2_match_rate':'NA'},cfg), 'ambiguous')
        self.assertEqual(classify({'n1_match_rate':'NA','n2_match_rate':.2,
                                   'd1_match_rate':.1,'d2_match_rate':'NA'},cfg), 'ambiguous')

    def test_summary_rows(self):
        classes = {'POP': [
            {'chromosome':'1','segment_id':'a','archaic_class':'neanderthal'},
            {'chromosome':'1','segment_id':'b','archaic_class':'ambiguous'}]}
        wide = {'POP': [
            {'chromosome':'1','segment_id':'a','length_bp':'1000000'},
            {'chromosome':'1','segment_id':'b','length_bp':'2000000'}]}
        result = classification_summary_rows(classes, wide, ['POP'])
        self.assertEqual(result[0]['n_segments'], 1)
        self.assertEqual(result[0]['genome_coverage_mb'], 1.0)
        calls = {'POP': [
            {'archaic_class':'denisovan','sample_id':'S1','length_bp':'1000000'},
            {'archaic_class':'denisovan','sample_id':'S1','length_bp':'2000000'},
            {'archaic_class':'neanderthal','sample_id':'S2','length_bp':'1000000'}]}
        samples = [{'sample_id':'S1','population':'POP','role':'target'},
                   {'sample_id':'S2','population':'POP','role':'target'}]
        summary = individual_summary_rows(calls, samples, ['POP'])
        deni = next(row for row in summary if row['archaic_class'] == 'denisovan')
        self.assertEqual(deni['per_individual_introgressed_mb'], 1.5)
        self.assertEqual(deni['n_introgressed_individuals'], 1)

    def test_collect_empty_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            p, q, out = [Path(tmp)/x for x in ('a.tsv','b.tsv','out.tsv.gz')]
            write(p, ['a'], [])
            write(q, ['a'], [{'a':1}])
            collect([p,q],out)
            self.assertEqual(read(out),[{'a':'1'}])


if __name__ == '__main__':
    unittest.main()
