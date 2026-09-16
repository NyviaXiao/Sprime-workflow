import csv, gzip, importlib.util, pathlib, tempfile, unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('affinity', ROOT/'workflow/scripts/build_affinity_tables.py')
affinity = importlib.util.module_from_spec(spec); spec.loader.exec_module(affinity)

class AffinityTests(unittest.TestCase):
    def test_projection(self):
        with tempfile.TemporaryDirectory() as d:
            d=pathlib.Path(d); wide=d/'wide.tsv'; cls=d/'class.tsv'; out=d/'a.tsv.gz'
            fields=['population','chromosome','segment_id','start','end','length_bp','altai_nean_matched','altai_nean_mismatch','altai_nean_callable','altai_nean_notcomp','altai_nean_match_rate']
            with wide.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerow(dict(zip(fields,['P','1','s1','10','20','11','3','1','4','0','.75'])))
            with cls.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=['chromosome','segment_id','archaic_class'],delimiter='\t'); w.writeheader(); w.writerow({'chromosome':'1','segment_id':'s1','archaic_class':'neanderthal'})
            affinity.build(wide, cls, 'P', 'altai_nean', out)
            with gzip.open(out,'rt',newline='') as f: row=next(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(row['match'],'3'); self.assertEqual(row['match_rate'],'.75'); self.assertEqual(row['archaic_class'],'neanderthal')

if __name__ == '__main__': unittest.main()
