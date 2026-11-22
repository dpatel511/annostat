"""Dataset-level hint when internal stops suggest a different NCBI translation table.

GFF3 coordinates are 1-based inclusive. Each constructed CDS is ATG, 20 GCT
codons, one internal test codon, 20 GCT codons and a terminal TAA (129 bp).
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from annostat.comparison import GenomeInput, run_comparison
from annostat.workflow import run_analysis

SPACER = "C" * 15


def build(root, internal_codons, extra_attributes=None):
    """Write one sequence carrying one CDS per internal codon ('GCT' = no internal stop)."""
    sequence, rows = "", []
    for index, codon in enumerate(internal_codons):
        start = len(sequence) + len(SPACER) + 1
        sequence += SPACER + "ATG" + "GCT" * 20 + codon + "GCT" * 20 + "TAA"
        extra = (extra_attributes or {}).get(index, "")
        rows.append(f"chr\ttest\tCDS\t{start}\t{start + 128}\t.\t+\t0\tID=cds{index}{extra}")
    fasta, gff = root / "genome.fna", root / "annotation.gff3"
    fasta.write_text(">chr\n" + sequence + SPACER + "\n", encoding="utf-8")
    gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return fasta, gff


class GeneticCodeCheckTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def analyze(self, codons, genetic_code=None, extra=None):
        fasta, gff = build(self.root, codons, extra)
        return run_analysis(fasta, gff, self.root / "out", "csv", genetic_code=genetic_code)

    def test_widespread_tga_stops_name_the_codon_and_every_reading_table(self):
        summary = self.analyze(["TGA"] * 12 + ["GCT"] * 8)
        check = summary["genetic_code_check"]
        self.assertEqual(check["models_with_internal_stops"], 12)
        self.assertEqual(check["cds_models"], 20)
        self.assertEqual(check["translation_tables"], [11])
        self.assertEqual(check["internal_stop_codons"], {"TGA": 12})
        self.assertEqual(check["alternative_tables"], [3, 4, 25])
        self.assertEqual(check["resolved_models"], 12)
        issues = summary["quality_control"]["issue_counts"]
        self.assertEqual(issues["possible_genetic_code_mismatch"], 1)
        self.assertEqual(issues["internal_stop_codon"], 12)
        report = (self.root / "out/report.html").read_text(encoding="utf-8")
        self.assertIn("Possible genetic-code mismatch", report)
        self.assertIn("TGA", report)

    def test_correct_table_has_no_hint(self):
        summary = self.analyze(["TGA"] * 12 + ["GCT"] * 8, genetic_code=4)
        self.assertIsNone(summary["genetic_code_check"])
        self.assertNotIn("possible_genetic_code_mismatch", summary["quality_control"]["issue_counts"])

    def test_thresholds_on_both_sides(self):
        cases = (
            (["TGA"] * 10 + ["GCT"] * 190, True),    # 10 models, 5.0% of CDS models
            (["TGA"] * 9 + ["GCT"] * 11, False),     # fewer than 10 models
            (["TGA"] * 10 + ["GCT"] * 191, False),   # under 5% of CDS models
            (["TGA"] * 9 + ["TAA"], True),           # 9 of 10 resolved by another table
            (["TGA"] * 8 + ["TAA"] * 2, False),      # 8 of 10 resolved
        )
        for codons, expected in cases:
            with self.subTest(codons=f"{codons.count('TGA')} TGA, {codons.count('TAA')} TAA, {len(codons)} CDS"):
                self.assertEqual(self.analyze(codons)["genetic_code_check"] is not None, expected)

    def test_stops_no_supported_table_reads_have_no_hint(self):
        summary = self.analyze(["TAA", "TAG"] * 6)
        self.assertIsNone(summary["genetic_code_check"])
        self.assertEqual(summary["quality_control"]["issue_counts"]["internal_stop_codon"], 12)

    def test_pseudogenes_and_declared_tables_follow_the_internal_stop_rule(self):
        pseudo = {index: ";pseudo=true" for index in range(3)}
        summary = self.analyze(["TGA"] * 15, extra=pseudo)
        self.assertEqual(summary["genetic_code_check"]["models_with_internal_stops"], 12)
        self.assertEqual(summary["quality_control"]["issue_counts"]["internal_stop_codon"], 12)
        declared = {index: ";transl_table=4" for index in range(6)}
        summary = self.analyze(["TGA"] * 24, extra=declared)
        check = summary["genetic_code_check"]
        self.assertEqual((check["models_with_internal_stops"], check["cds_models"]), (18, 24))
        self.assertEqual(check["translation_tables"], [11])

    def test_aggregate_comparison_counts_the_same_finding(self):
        fasta, gff = build(self.root, ["TGA"] * 12 + ["GCT"] * 8)
        result = run_comparison([GenomeInput("A", fasta, gff), GenomeInput("B", fasta, gff)], self.root / "cmp")
        for profile in result["datasets"]:
            self.assertEqual(profile["qc_issue_counts"]["possible_genetic_code_mismatch"], 1)

    def test_terminal_summary_names_the_codon_and_quiet_suppresses_it(self):
        fasta, gff = build(self.root, ["TGA"] * 12 + ["GCT"] * 8)
        command = [sys.executable, "-m", "annostat", "analyze", "-f", str(fasta), "-g", str(gff)]
        result = subprocess.run([*command, "-o", str(self.root / "cli")], capture_output=True, text=True, check=True)
        self.assertIn("Genetic code check", result.stdout)
        self.assertIn("12/20 CDS models stop internally at TGA; tables 3, 4, 25 translate it", result.stdout)
        quiet = subprocess.run([*command, "-o", str(self.root / "quiet"), "-q"], capture_output=True, text=True, check=True)
        self.assertEqual((quiet.stdout, quiet.stderr), ("", ""))
        written = json.loads((self.root / "cli/summary.json").read_text(encoding="utf-8"))
        self.assertEqual(written["genetic_code_check"]["alternative_tables"], [3, 4, 25])


if __name__ == "__main__":
    unittest.main()
