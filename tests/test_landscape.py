import importlib.util
import pathlib
import shutil
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).parents[1]


class LandscapeSmokeTests(unittest.TestCase):
    def test_scripts_normalize_chromosomes_and_use_named_outputs(self):
        intro = (ROOT/'workflow/scripts/plot_landscape.R').read_text()
        affinity = (ROOT/'workflow/scripts/plot_affinity_landscape.R').read_text()
        for script in (intro, affinity):
            self.assertIn('normalize_chr <- function(x)', script)
            self.assertIn('ifelse(grepl("^chr", x), x, paste0("chr", x))', script)
            self.assertIn('genome_build', script)
        self.assertIn('output[["introgression"]]', intro)
        self.assertIn('plot.type=6', intro)
        self.assertIn('data.panel="ideogram"', intro)
        self.assertNotIn('data.panel=1', intro)
        self.assertIn('y0=0, y1=1', intro)
        self.assertNotIn('y0=0.18', intro)
        self.assertIn('kpRect(kp, chr=chromosomes', intro)
        self.assertIn('legend("bottom"', intro)
        self.assertNotIn('kpAddLabels', affinity)
        self.assertIn('plot.type=6', affinity)
        self.assertIn('data.panel="ideogram"', affinity)
        self.assertIn('lower <- (n - i) / n', affinity)
        self.assertIn('upper <- (n - i + 1) / n', affinity)
        self.assertIn('add_colorbar <- function', affinity)
        self.assertIn('rasterImage', affinity)

    def test_r_landscape_smoke_when_runtime_is_available(self):
        rscript = shutil.which('Rscript')
        if not rscript:
            self.skipTest('Rscript is not installed')
        package = subprocess.run([rscript, '-e', 'library(karyoploteR)'], capture_output=True)
        if package.returncode:
            self.skipTest('karyoploteR is not installed')
        with tempfile.TemporaryDirectory() as tmp:
            tmp = pathlib.Path(tmp)
            classification = tmp/'classification.tsv'
            classification.write_text('chromosome\tstart\tend\tarchaic_class\n1\t10\t20\tneanderthal\nchr1\t30\t40\tdenisovan\n')
            affinity = tmp/'n1.tsv'; affinity.write_text('chromosome\tstart\tend\tarchaic_class\tmatch_rate\n1\t10\t20\tneanderthal\t0.8\n')
            intro_png, nean_png, deni_png = tmp/'intro.png', tmp/'nean.png', tmp/'deni.png'
            wrapper = tmp/'smoke.R'
            wrapper.write_text(f'''setClass("Snakemake", slots=c(input="list", output="list", params="list", wildcards="list"))
snakemake <- new("Snakemake", input=list(classification="{classification.as_posix()}", nean_aff=list("{affinity.as_posix()}"), deni_aff=list("{affinity.as_posix()}")), output=list(introgression="{intro_png.as_posix()}", neanderthal="{nean_png.as_posix()}", denisovan="{deni_png.as_posix()}"), params=list(chromosomes=list("1", "chr1"), genome_build="GRCh37", neanderthal_refs=list("n1"), denisovan_refs=list("d1")), wildcards=list(population="P"))
source("{(ROOT/'workflow/scripts/plot_landscape.R').as_posix()}")
source("{(ROOT/'workflow/scripts/plot_affinity_landscape.R').as_posix()}")
''')
            completed = subprocess.run([rscript, str(wrapper)], capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            for path in (intro_png, nean_png, deni_png):
                self.assertGreater(path.stat().st_size, 0)


if __name__ == '__main__': unittest.main()
