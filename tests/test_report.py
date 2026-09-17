import base64
import csv
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
        with open(path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
            writer.writeheader(); writer.writerows(rows)

    def test_research_report_html_and_pdf(self):
        try:
            import matplotlib, weasyprint
        except ImportError:
            self.skipTest('matplotlib and WeasyPrint are required for report integration')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp); report_dir = tmp/'report'; report_dir.mkdir()
            run = tmp/'run.json'; run.write_text(json.dumps({'run_time':'now','genome_build':'GRCh37','populations':['P'],'chromosomes':['1'],'archaic_references':['n1','d1'],'enabled_modules':{'report':True,'contour':True},'git_commit':'abc'}))
            software = tmp/'software.tsv'; self.write_tsv(software, ['software','version'], [{'software':'python','version':'x'}])
            classes = tmp/'classification.tsv'; self.write_tsv(classes, ['population','archaic_class','n_segments','genome_coverage_mb'], [{'population':'P','archaic_class':'neanderthal','n_segments':10,'genome_coverage_mb':1.2},{'population':'P','archaic_class':'denisovan','n_segments':5,'genome_coverage_mb':.5},{'population':'P','archaic_class':'ambiguous','n_segments':2,'genome_coverage_mb':.2}])
            individual = tmp/'individual.tsv'; self.write_tsv(individual, ['population','archaic_class','per_individual_introgressed_mb','n_introgressed_individuals'], [{'population':'P','archaic_class':'neanderthal','per_individual_introgressed_mb':1,'n_introgressed_individuals':1}])
            gmm = tmp/'gmm.tsv'; self.write_tsv(gmm, ['population','reference','selected_components'], [{'population':'P','reference':'d1','selected_components':'1'},{'population':'P','reference':'d2','selected_components':'2'}])
            adaptive = tmp/'top2.tsv'; self.write_tsv(adaptive, ['population','adaptive_group','adaptive_rank','chromosome','segment_id','candidate_start','candidate_end','core_mean_AF','max_AF','core_variant_count','neanderthal_pass','denisovan_pass'], [{'population':'P','adaptive_group':'neanderthal','adaptive_rank':1,'chromosome':'1','segment_id':'n','candidate_start':1,'candidate_end':2,'core_mean_AF':.8,'max_AF':.8,'core_variant_count':1,'neanderthal_pass':1,'denisovan_pass':0},{'population':'P','adaptive_group':'denisovan','adaptive_rank':1,'chromosome':'1','segment_id':'d','candidate_start':3,'candidate_end':4,'core_mean_AF':.7,'max_AF':.7,'core_variant_count':1,'neanderthal_pass':0,'denisovan_pass':1}])
            shared = tmp/'shared.tsv'; self.write_tsv(shared, ['adaptive_group','chromosome','start','end','populations','n_populations'], [])
            one, two, contour = tmp/'P.d1.png', tmp/'P.d2.png', tmp/'P.n1__d1.png'
            for path in (one, two, contour): path.write_bytes(PNG)
            html_path, pdf_path = report_dir/'report.html', report_dir/'report.pdf'
            report.build_report(html_path, pdf_path, run, software, classes, [], [], contour, individual, gmm, [one, two], adaptive, shared, report_dir/'classification_summary.png', report_dir/'individual_summary.png')
            page = html_path.read_text(encoding='utf-8')
            for heading in ('Executive Summary', 'A. Run Overview', 'B. Classification Results', 'C. Introgression Landscape', 'D1. Neanderthal Affinity', 'D2. Denisovan Affinity', 'F. Individual-level Introgression', 'G. GMM Analysis', 'H. Adaptive Introgression Candidates', 'I. Reproducibility & Provenance'):
                self.assertIn(heading, page)
            self.assertIn('10 Neanderthal-classified', page)
            self.assertIn('P.d2.png', page); self.assertNotIn('P.d1.png', page)
            self.assertNotIn('resolved_config', page)
            self.assertGreater(pdf_path.stat().st_size, 0)

    def test_disabled_optional_inputs_do_not_require_files(self):
        try:
            import matplotlib, weasyprint
        except ImportError:
            self.skipTest('matplotlib and WeasyPrint are required for report integration')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp); run = tmp/'run.json'; run.write_text(json.dumps({'run_time':'now','genome_build':'GRCh37','populations':['P'],'chromosomes':['1'],'archaic_references':[],'enabled_modules':{},'git_commit':'abc'}))
            software = tmp/'software.tsv'; self.write_tsv(software, ['software','version'], [])
            classes = tmp/'classification.tsv'; self.write_tsv(classes, ['population','archaic_class','n_segments','genome_coverage_mb'], [{'population':'P','archaic_class':'neanderthal','n_segments':0,'genome_coverage_mb':0}])
            html_path, pdf_path = tmp/'report'/'report.html', tmp/'report'/'report.pdf'
            report.build_report(html_path, pdf_path, run, software, classes, [], [], None, None, None, [], None, None, html_path.parent/'classification_summary.png')
            self.assertIn('Not enabled for this run.', html_path.read_text())
            self.assertTrue(pdf_path.exists())


if __name__ == '__main__': unittest.main()
