import base64
import csv
import gzip
import importlib.util
import json
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('report', ROOT/'workflow/scripts/build_report.py')
report = importlib.util.module_from_spec(spec); spec.loader.exec_module(report)
PNG = base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVQIHWP4z8DwHwAFgAI/ScL0rQAAAABJRU5ErkJggg==')


class ReportTests(unittest.TestCase):
    def write_tsv(self, path, fields, rows):
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'wt', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader(); writer.writerows(rows)

    def config(self):
        return {'sprime_minscore':150000, 'summary_min_callable':10,
                'archaic_references':[{'id':'n1','tag':'Altai','group':'neanderthal'}, {'id':'d1','tag':'Denisovan3','group':'denisovan'}],
                'classification':{'neanderthal_references':['n1'],'denisovan_references':['d1'],'neanderthal_gt':.6,'denisovan_lt_for_neanderthal':.4,'denisovan_gt':.3,'neanderthal_lt_for_denisovan':.3},
                'gmm':{'enabled':True,'score_gt':150000,'min_callable':30,'min_segments':10,'method':'legacy_chi_square','alpha':.05,'target_references':['d1']},
                'individual':{'enabled':True,'min_markers':2},
                'adaptive':{'enabled':True,'min_allele_frequency':.3,'max_frequency_drop':.2,'min_callable_sites':10,'neanderthal_match_threshold':.5,'denisovan_match_threshold':.4},
                'report':{'contour_population':'P','contour_pair':['n1','d1']}}

    def test_narratives_distinguish_screening_from_top2(self):
        segment_rows = [{'neanderthal_pass':'1','denisovan_pass':'0'}] * 4 + [{'neanderthal_pass':'0','denisovan_pass':'1'}] * 3
        top2_rows = [{'population':'P','adaptive_group':'neanderthal'}] * 2 + [{'population':'P','adaptive_group':'denisovan'}] * 2
        text = report._adaptive_narrative(segment_rows, top2_rows)
        self.assertIn('4 Neanderthal-associated and 3 Denisovan-associated passing segments', text)
        self.assertNotIn('candidate rows were retained', text)
        individual = report._individual_narrative([{'population':'P','archaic_class':'neanderthal','per_individual_introgressed_mb':'4.21','n_introgressed_individuals':'103'}])
        self.assertIn('4.21 Mb', individual); self.assertIn('103 individuals', individual)

    def test_research_report_html_and_pdf(self):
        try:
            import matplotlib, weasyprint
        except ImportError:
            self.skipTest('matplotlib and WeasyPrint are required for report integration')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp); report_dir = tmp/'report'; report_dir.mkdir()
            config = tmp/'config.json'; config.write_text(json.dumps(self.config()))
            run = tmp/'run.json'; run.write_text(json.dumps({'run_time':'now','genome_build':'GRCh37','populations':['P'],'chromosomes':['1'],'archaic_references':['n1','d1'],'enabled_modules':{'report':True},'git_commit':'abc','git_dirty':False,'git_diff_sha256':'digest'}))
            software = tmp/'software.tsv'; self.write_tsv(software, ['software','version'], [{'software':'python','version':'x'}])
            classes = tmp/'classification.tsv'; self.write_tsv(classes, ['population','archaic_class','n_segments','genome_coverage_mb'], [{'population':'P','archaic_class':'neanderthal','n_segments':10,'genome_coverage_mb':1.2},{'population':'P','archaic_class':'denisovan','n_segments':5,'genome_coverage_mb':.5},{'population':'P','archaic_class':'ambiguous','n_segments':2,'genome_coverage_mb':.2}])
            individual = tmp/'individual.tsv'; self.write_tsv(individual, ['population','archaic_class','per_individual_introgressed_mb','n_introgressed_individuals'], [{'population':'P','archaic_class':'neanderthal','per_individual_introgressed_mb':4.21,'n_introgressed_individuals':103}])
            gmm = tmp/'gmm.tsv'; self.write_tsv(gmm, ['population','reference','n_segments','status','method','selected_components','adjusted_p_1_vs_2','second_comparison','adjusted_p_second','converged','loglik_1','aic_1','bic_3'], [{'population':'P','reference':'d1','n_segments':10,'status':'ok','method':'legacy','selected_components':'2','adjusted_p_1_vs_2':'.01','second_comparison':'2_vs_3','adjusted_p_second':'.2','converged':'True','loglik_1':'bad','aic_1':'bad','bic_3':'bad'}])
            segments = tmp/'segments.tsv'; self.write_tsv(segments, ['population','neanderthal_pass','denisovan_pass'], [{'population':'P','neanderthal_pass':1,'denisovan_pass':0}]*4 + [{'population':'P','neanderthal_pass':0,'denisovan_pass':1}]*3)
            top2 = tmp/'top2.tsv'; top_fields = ['population','adaptive_group','adaptive_rank','chromosome','segment_id','candidate_start','candidate_end','core_mean_AF','max_AF','core_variant_count','neanderthal_pass','denisovan_pass']; self.write_tsv(top2, top_fields, [{'population':'P','adaptive_group':'neanderthal','adaptive_rank':1,'chromosome':'1','segment_id':'n','candidate_start':1,'candidate_end':2,'core_mean_AF':.8,'max_AF':.8,'core_variant_count':1,'neanderthal_pass':1,'denisovan_pass':0},{'population':'P','adaptive_group':'denisovan','adaptive_rank':1,'chromosome':'1','segment_id':'d','candidate_start':3,'candidate_end':4,'core_mean_AF':.7,'max_AF':.7,'core_variant_count':1,'neanderthal_pass':0,'denisovan_pass':1}])
            shared = tmp/'shared.tsv'; self.write_tsv(shared, ['adaptive_group','chromosome','start','end','populations','n_populations'], [])
            affinity = tmp/'P'/'n1.tsv.gz'; affinity.parent.mkdir(); self.write_tsv(affinity, ['population','archaic_class','match_rate'], [{'population':'P','archaic_class':'neanderthal','match_rate':.6},{'population':'P','archaic_class':'neanderthal','match_rate':.8}])
            one, contour = tmp/'P.d1.png', tmp/'P.n1__d1.png'
            for path in (one, contour): path.write_bytes(PNG)
            html_path, pdf_path = report_dir/'report.html', report_dir/'report.pdf'
            report.build_report(html_path, pdf_path, run, software, config, classes, [], [], [affinity], contour, individual, gmm, [one], segments, top2, shared, report_dir/'classification_summary.png', report_dir/'individual_summary.png')
            page = html_path.read_text(encoding='utf-8')
            for heading in ('Executive Summary', 'A. Run Overview', 'A2. Analysis Settings', 'B. Classification Results', 'D1. Neanderthal Affinity', 'D2. Denisovan Affinity', 'I. Reproducibility & Provenance'):
                self.assertIn(heading, page)
            self.assertIn('10 Neanderthal-classified', page); self.assertIn('4.21 Mb', page)
            self.assertIn('4 Neanderthal-associated and 3 Denisovan-associated passing segments', page)
            self.assertIn('Altai (n=2; median=0.70; mean=0.70)', page)
            self.assertIn('selected_components', page); self.assertIn('adjusted_p_1_vs_2', page); self.assertIn('adjusted_p_second', page)
            self.assertNotIn('loglik_1', page); self.assertNotIn('aic_1', page); self.assertNotIn('bic_3', page)
            self.assertGreater(pdf_path.stat().st_size, 0)


if __name__ == '__main__': unittest.main()
