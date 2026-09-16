import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('provenance', ROOT/'workflow/scripts/build_provenance.py')
provenance = importlib.util.module_from_spec(spec); spec.loader.exec_module(provenance)

class ProvenanceTests(unittest.TestCase):
    def test_declared_inputs_expand_every_chromosome(self):
        cfg={'samples':'samples.tsv','genetic_map':'map','sprime_jar':'jar','map_arch':'arch','vcf':'modern.{chrom}.vcf.gz','chromosomes':['1','2'],'archaic_references':[{'id':'n1','vcf':'n1.{chrom}.vcf.gz','mask':'n1.{chrom}.bed.gz'}]}
        rows=dict(provenance._declared_inputs(cfg))
        self.assertEqual(rows['vcf_chr1'],'modern.1.vcf.gz'); self.assertEqual(rows['vcf_chr2'],'modern.2.vcf.gz')
        self.assertEqual(rows['n1_mask_chr2'],'n1.2.bed.gz')

if __name__ == '__main__': unittest.main()
