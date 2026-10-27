"""Evidence consistency and safe offline rendering of the complete CDS report."""
import csv
import json
import re
import tempfile
import unittest
from pathlib import Path

from annostat.comparison import GenomeInput, run_comparison


class GeneReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        fasta = self.root / "genome.fna"
        # 1-based inclusive: long CDS 1..21 = MAMEQQ*, short 7..21 = MEQQ*.
        fasta.write_text(">chr\nATGGCTATGGAACAACAATAA" + "ATG" * 400 + "NNNATGAAATAAGTGAAACAATAA\n", encoding="utf-8")
        self.datasets = []
        for label, start in (("A", 1), ("B", 7)):
            gff = self.root / f"{label}.gff3"
            rows = [f"chr\ttest\tCDS\t{start}\t21\t.\t+\t0\tID=gyrB_{label};gene=gyrB;product=DNA gyrase"]
            # More than the old 100-row preview limit, including hostile metadata.
            rows += [f"chr\ttest\tCDS\t{30+i*9}\t{38+i*9}\t.\t+\t0\tID=cds_{i};product=%3C/script%3E%3Cimg src=x onerror=alert(1)%3E" for i in range(110)]
            rows.append(f"chr\ttest\tCDS\t{1100 if label == 'A' else 1150}\t{1120 if label == 'A' else 1170}\t.\t+\t0\tID=only_{label}")
            rows.append(f"chr\ttest\tCDS\t1222\t1233\t.\t+\t0\tID=limited_{label}")
            rows.append(f"chr\ttest\tCDS\t1234\t1245\t.\t+\t0\tID=initiation_{label}" + (";partial=true" if label == "B" else ""))
            gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
            self.datasets.append(GenomeInput(label, fasta, gff))

    def test_all_outputs_share_measured_consequences_and_complete_safe_payload(self):
        output = self.root / "out"
        saved = run_comparison(self.datasets, output, gene_level=True)
        gene = saved["gene_level_comparison"]
        row = next(r for r in gene["relationships"] if r["id_a"] == "gyrB_A")
        self.assertEqual(row["protein_consequence"], "terminal_only")
        self.assertEqual((row["protein_length_delta_aa"], row["n_terminal_delta_aa"]), (-2, -2))
        self.assertEqual((row["n_terminal_residues_a"], row["identical_aa"]), ("MA", 4))
        self.assertEqual((row["start_delta_bp"], row["start_codon_a"], row["start_codon_b"]), (6, "ATG", "ATG"))
        self.assertEqual(sum(gene["consequence_counts"].values()), gene["status_counts"]["matched"])
        self.assertEqual(gene["unmatched_cds_counts"], {"A": 1, "B": 1})
        html = (output / "comparison.html").read_text(encoding="utf-8")
        payload = json.loads(re.search(r'<script id="gene-data" type="application/json">(.*?)</script>', html, re.S).group(1))
        self.assertEqual(payload["relationships"], gene["relationships"])
        self.assertGreater(len(payload["relationships"]), 100)
        self.assertNotIn("</script><img", html)
        self.assertNotRegex(html, r'<script[^>]+src=|<link[^>]+href="https?://')
        self.assertIn("DataTables 3.1.4", html)
        self.assertIn("comparison.txt", html)
        self.assertIn("1-based inclusive", html)
        with (output / "tables/gene_comparisons.tsv").open(newline="", encoding="utf-8") as handle:
            tsv = next(r for r in csv.DictReader(handle, delimiter="\t") if r["id_a"] == "gyrB_A")
        self.assertEqual(tsv["n_terminal_residues_a"], row["n_terminal_residues_a"])
        self.assertEqual(int(tsv["protein_length_delta_aa"]), -2)
        self.assertEqual(json.loads((output / "comparison.json").read_text(encoding="utf-8")), saved)
        text = (output / "comparison.txt").read_text(encoding="utf-8")
        self.assertIn("Resolved one-to-one CDS pairs: 113", text)
        self.assertIn("Terminal-only difference: 1", text)
        self.assertIn("gyrB_A", text)
        self.assertIn("MA", text)
        self.assertIn("4 aa", text)
        # A global length ranking hid unknown translations and unmatched calls.
        for title, count in (("Terminal boundary differences", 1), ("Translation-context differences", 1),
                             ("Limited protein evidence", 1), ("Unresolved correspondence", 2)):
            caption = f"{title}: {count} relationship row{'s' if count != 1 else ''}"
            self.assertIn(caption, text)
        # Displayed percentages use accepted pairs, not all relationship rows.
        self.assertIn('110 of 113 accepted pairs (97.3%)', html)
        self.assertIn('2 of 113 accepted pairs (1.8%)', html)
        self.assertIn('1 of 113 accepted pairs (0.9%)', html)
        self.assertIn('113 / 114 CDS models (99.1%)', html)
        self.assertIn('Exact CDS structure', html)
        self.assertIn('Alternative start boundary', html)
        self.assertIn('2 unresolved relationship rows', html)
        self.assertIn('B shorter: <b>1</b> · B longer: <b>0</b> · Same length: <b>1</b>', html)
        self.assertIn('median absolute Δ: <strong>1 aa</strong>', html)
        self.assertIn('signed range: -2 to +0 aa', html)
        self.assertIn('2 × 112 shared models / 228 CDS models', html)
        self.assertIn('id="gene-overview"', html)
        self.assertIn('id="gene-differences"', html)
        self.assertIn('id="gene-evidence"', html)
        for identifier in ("limited_A", "initiation_A", "only_A", "only_B"):
            self.assertIn(identifier, text)
        self.assertNotIn("None aa", text)
        self.assertIn("unknown/internal-stop residues", text)
        self.assertIn("B-A +0 aa", text)
        svg = (output / "plots/gene_consequences.svg").read_text(encoding="utf-8")
        self.assertIn("113 resolved", svg)
        self.assertIn("Terminal-only", svg)
        run_comparison(self.datasets, output)
        self.assertFalse((output / "comparison.txt").exists())
        self.assertFalse((output / "plots/gene_consequences.svg").exists())

    def test_empty_gene_report_keeps_zero_and_unavailable_distinct(self):
        for dataset in self.datasets:
            dataset.gff.write_text("##gff-version 3\n", encoding="utf-8")
        output = self.root / "empty"
        result = run_comparison(self.datasets, output, gene_level=True)["gene_level_comparison"]
        self.assertEqual(result["consequence_counts"], {})
        self.assertIn("Resolved one-to-one CDS pairs: 0", (output / "comparison.txt").read_text(encoding="utf-8"))
        self.assertIn("n/a", (output / "comparison.txt").read_text(encoding="utf-8"))
        html = (output / "comparison.html").read_text(encoding="utf-8")
        self.assertIn('No one-to-one CDS pairs were resolved', html)
        self.assertNotIn('0.0% of accepted pairs', html)
        self.assertIn('Not available', html)


if __name__ == "__main__":
    unittest.main()
