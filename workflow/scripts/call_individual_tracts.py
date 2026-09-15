"""Call individual haplotype runs within population SPrime segments.

This replaces the former shell/AWK/paste/selectGene2 chain. VCF and SPrime
records are joined by CHROM, POS, REF and ALT before haplotypes are examined.
"""
from collections import defaultdict
from core import CALLS, read, write, score_rows, runs


def individuals(s, c):
    """Write individual carried-marker runs for one population/chromosome."""
    import pysam
    # Exclude outgroup samples at input rather than removing them from a large
    # result file after calling.
    targets = [r['sample_id'] for r in read(s.input.samples)
               if r['role'] == 'target'
               and r['population'] == s.wildcards.population]
    labels = {(r['chromosome'], r['segment_id']): r['archaic_class'] for r in read(s.input.classes)}
    # Group markers by SPrime segment. Classification is segment-level, while
    # carrying status is evaluated separately for every sample and haplotype.
    segments = defaultdict(list)
    for r in score_rows(s.input.score):
        segments[r['CHROM'], r['SEGMENT']].append(r)
    with pysam.VariantFile(s.input.vcf) as v:
        def emit():
            for (chrom, seg), rows in segments.items():
                rows.sort(key=lambda r: int(r['POS']))
                # Alleles are part of the key: position-only joins can silently
                # match incompatible variants after normalization changes.
                expected = {}
                for r in rows:
                    key = (int(r['POS']), r['REF'], r['ALT'])
                    if key in expected:
                        raise ValueError(f'Duplicate SPrime marker {key}')
                    expected[key] = r
                # Fetch only the span occupied by this segment, then retain the
                # exact SPrime markers. CSI is produced by the subset rule.
                genotypes = {}
                for rec in v.fetch(chrom, int(rows[0]['POS']) - 1, int(rows[-1]['POS'])):
                    key = (rec.pos, rec.ref, ','.join(rec.alts or []))
                    if key not in expected:
                        continue
                    if key in genotypes:
                        raise ValueError(f'Duplicate VCF record {key}')
                    calls = {}
                    for sample in targets:
                        call = rec.samples[sample]
                        gt = call.get('GT')
                        if gt is None or len(gt) != 2:
                            raise ValueError(f'Diploid GT required: {sample} {chrom}:{rec.pos}')
                        # Missing alleles are allowed and will break a marker
                        # run. Any observed diploid call must be phased.
                        if not call.phased and any(a is not None for a in gt):
                            raise ValueError(f'Phased GT required: {sample} {chrom}:{rec.pos}')
                        calls[sample] = gt
                    genotypes[key] = calls
                if set(expected) != set(genotypes):
                    raise ValueError(f'SPrime markers missing from VCF: {chrom} {seg}')
                # ALLELE is the SPrime allele index (0=REF, 1=ALT). A haplotype
                # carries the candidate marker when its GT allele equals it.
                for sample in targets:
                    for hap in (0, 1):
                        markers = []
                        for key, row in expected.items():
                            allele = int(row['ALLELE'])
                            if allele not in (0, 1):
                                raise ValueError('Only biallelic SPrime alleles supported')
                            markers.append((key[0], genotypes[key][sample][hap] == allele))
                        for start, end, count in runs(markers, c['individual']['min_markers']):
                            yield dict(zip(CALLS, [s.wildcards.population, sample, hap + 1, chrom, start, end,
                                                   end-start+1, count, seg, labels[(chrom, seg)]]))
        write(s.output[0], CALLS, emit())
