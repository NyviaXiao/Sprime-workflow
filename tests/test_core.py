import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from core import classify, collect, gmm_pass, read, runs, summarize, summary_fields, write


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
        self.assertFalse(gmm_pass({**row,'d2_match_rate':'NA'},cfg))

    def test_classification_and_missing(self):
        cfg = dict(neanderthal_reference='n', denisovan_reference='d', neanderthal_gt=.6,
                   denisovan_lt_for_neanderthal=.4, denisovan_gt=.3, neanderthal_lt_for_denisovan=.3)
        for n,d,label in [(.8,.1,'neanderthal'),(.1,.8,'denisovan'),(.5,.5,'ambiguous'),
                          (.1,.1,'ambiguous'),('NA',.5,'ambiguous')]:
            self.assertEqual(classify({'n_match_rate':n,'d_match_rate':d},cfg),label)

    def test_collect_empty_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            p, q, out = [Path(tmp)/x for x in ('a.tsv','b.tsv','out.tsv.gz')]
            write(p, ['a'], [])
            write(q, ['a'], [{'a':1}])
            collect([p,q],out)
            self.assertEqual(read(out),[{'a':'1'}])


if __name__ == '__main__':
    unittest.main()
