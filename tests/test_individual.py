import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow/scripts'))
from call_individual_tracts import individuals
from core import read, write


class IndividualTests(unittest.TestCase):
    def test_vcf_alignment_haplotypes_and_phase(self):
        class Call(dict):
            phased = True

        class VCF:
            def __enter__(self):
                return self
            def __exit__(self, *args):
                pass
            def fetch(self, chrom, start, end):
                return iter(records)

        records = [SimpleNamespace(pos=pos, ref='A', alts=('G',),
                    samples={'sample.with.dot': Call(GT=gt)})
                   for pos,gt in [(10,(1,0)),(20,(1,0)),(30,(None,None)),(40,(0,1)),(50,(0,1))]]
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            write(p/'samples.tsv', ['sample_id','population','role'], [
                {'sample_id':'sample.with.dot','population':'POP','role':'target'},
                {'sample_id':'OUT','population':'YRI','role':'outgroup'}])
            write(p/'classes.tsv', ['chromosome','segment_id','archaic_class'],
                  [{'chromosome':'1','segment_id':'7','archaic_class':'denisovan'}])
            fields = ['CHROM','POS','ID','REF','ALT','SEGMENT','ALLELE','SCORE']
            write(p/'score',fields,[dict(zip(fields,['1',pos,'.','A','G','7','1',200000])) for pos in (10,20,30,40,50)])
            s = SimpleNamespace(input=SimpleNamespace(samples=p/'samples.tsv', classes=p/'classes.tsv',
                                score=p/'score',vcf='mock.vcf'), output=[p/'calls.tsv'],
                                wildcards=SimpleNamespace(population='POP'))
            with patch.dict(sys.modules, {'pysam':SimpleNamespace(VariantFile=lambda _:VCF())}):
                individuals(s, {'individual':{'min_markers':2}})
                calls = read(p/'calls.tsv')
                self.assertEqual([(r['start'],r['end'],r['haplotype']) for r in calls], [('10','20','1'),('40','50','2')])
                self.assertTrue(all(r['sample_id']=='sample.with.dot' and r['archaic_class']=='denisovan' for r in calls))
                records[0].samples['sample.with.dot'].phased = False
                with self.assertRaisesRegex(ValueError,'Phased GT required'):
                    individuals(s, {'individual':{'min_markers':2}})


if __name__ == '__main__':
    unittest.main()
