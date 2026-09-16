import csv
import importlib.util
import json
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('report', ROOT/'workflow/scripts/build_report.py')
report = importlib.util.module_from_spec(spec); spec.loader.exec_module(report)

class ReportTests(unittest.TestCase):
    def write_tsv(self, path, fields, rows):
        with open(path, 'w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t'); writer.writeheader(); writer.writerows(rows)

    def test_real_html_sections_and_selected_gmm_figures(self):
        try:
            import matplotlib  # report barplots are part of the integration.
        except ImportError:
            self.skipTest('matplotlib is not installed')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp); report_dir = tmp/'report'; report_dir.mkdir()
            run = tmp/'run.json'; run.write_text(json.dumps({'run_time':'now','genome_build':'GRCh37','populations':['P'],'chromosomes':['1'],'archaic_references':['n1','d1'],'enabled_modules':{'report':True},'git_commit':'abc'}))
            software=tmp/'software.tsv'; self.write_tsv(software,['software','version'],[{'software':'python','version':'x'}])
            classes=tmp/'classification.tsv'; self.write_tsv(classes,['population','archaic_class','n_segments','genome_coverage_mb'],[{'population':'P','archaic_class':'neanderthal','n_segments':1,'genome_coverage_mb':1},{'population':'P','archaic_class':'denisovan','n_segments':0,'genome_coverage_mb':0},{'population':'P','archaic_class':'ambiguous','n_segments':0,'genome_coverage_mb':0}])
            individual=tmp/'individual.tsv'; self.write_tsv(individual,['population','archaic_class','per_individual_introgressed_mb','n_introgressed_individuals'],[{'population':'P','archaic_class':'neanderthal','per_individual_introgressed_mb':1,'n_introgressed_individuals':1}])
            gmm=tmp/'gmm.tsv'; self.write_tsv(gmm,['population','reference','selected_components'],[{'population':'P','reference':'d1','selected_components':'1'},{'population':'P','reference':'d2','selected_components':'2'}])
            adaptive=tmp/'top2.tsv'; self.write_tsv(adaptive,['population','chromosome','candidate_start','candidate_end','core_mean_AF'],[{'population':'P','chromosome':'1','candidate_start':1,'candidate_end':2,'core_mean_AF':.8}])
            one, two, contour = tmp/'P.d1.png', tmp/'P.d2.png', tmp/'P.n1__d1.png'
            for path in (one,two,contour): path.write_bytes(b'png')
            output=report_dir/'report.html'
            report.build_report(output, run, software, classes, [], [], contour, individual, gmm, [one,two], adaptive)
            page=output.read_text()
            positions=[page.index(f'<h2>{letter}. {title}</h2>') for letter,title in zip('ABCDEFGH',report.SECTION_ORDER)]
            self.assertEqual(positions, sorted(positions)); self.assertIn('P.d2.png', page); self.assertNotIn('P.d1.png', page)
            self.assertIn('candidate_start', page); self.assertNotIn('resolved_config', page)

    def test_disabled_optional_inputs_do_not_require_files(self):
        try:
            import matplotlib
        except ImportError:
            self.skipTest('matplotlib is not installed')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp); run = tmp/'run.json'; run.write_text(json.dumps({'run_time':'now','genome_build':'GRCh37','populations':['P'],'chromosomes':['1'],'archaic_references':[],'enabled_modules':{},'git_commit':'abc'}))
            software=tmp/'software.tsv'; self.write_tsv(software,['software','version'],[])
            classes=tmp/'classification.tsv'; self.write_tsv(classes,['population','archaic_class','n_segments','genome_coverage_mb'],[{'population':'P','archaic_class':'neanderthal','n_segments':0,'genome_coverage_mb':0}])
            output=tmp/'report'/'report.html'
            report.build_report(output, run, software, classes, [], [], None, None, None, [], None)
            self.assertIn('Not enabled.', output.read_text())
            self.assertFalse((output.parent/'individual_summary.png').exists())

if __name__ == '__main__': unittest.main()
