import importlib.util, pathlib, unittest

ROOT = pathlib.Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('report', ROOT/'workflow/scripts/build_report.py')
report = importlib.util.module_from_spec(spec); spec.loader.exec_module(report)

class ReportTests(unittest.TestCase):
    def test_fixed_section_order(self):
        self.assertEqual(report.SECTION_ORDER, ['Run overview','Classification summary','Introgression landscape','Landscape + affinity','Neanderthal × Denisovan contour plot','Individual summary','GMM summary','Adaptive candidate summary'])

if __name__ == '__main__': unittest.main()
