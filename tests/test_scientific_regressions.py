"""Scientific audit regressions; genomic coordinates are 1-based inclusive."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from annostat.cohort import build_cohort, write_cohort
from annostat.comparison import GenomeInput, _annotation_profile, _pairwise_rows, run_comparison
from annostat.workflow import run_analysis
from annostat.gene_comparison import cds_model
from annostat.models import Feature
from annostat.sequences import iter_cds_sequences


class ScientificRegressionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.fasta = self.root / "genome.fna"
        self.fasta.write_text(">chr\nATGGTTGTGATGCAATAA\n", encoding="utf-8")

    def dataset(self, label, rows):
        gff = self.root / f"{label}.gff3"
        gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return GenomeInput(label, self.fasta, gff)

    def cds(self, identifier, start=1, end=18, strand="+", phase=0, attributes=""):
        return f"chr\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID={identifier}{attributes}"

    def compare_row(self, a, b):
        return run_comparison([a, b], self.root / "comparison", gene_level=True)["gene_level_comparison"]["relationships"][0]

    def test_coordinate_mapping_explains_initiator_recoding_without_moving_gaps(self):
        a = self.dataset("A", [self.cds("a")])
        b = self.dataset("B", [self.cds("b", 7)])
        row = self.compare_row(a, b)
        self.assertEqual(row["protein_consequence"], "initiator_recoding")
        self.assertEqual((row["alignment_a"], row["alignment_b"]), ("MVVMQ", "--MMQ"))
        self.assertEqual((row["identical_aa"], row["substitution_count"], row["gap_aa"]), (2, 1, 2))
        self.assertEqual((row["initiator_mismatch_count"], row["internal_mismatch_count"], row["shared_codon_count"]), (1, 0, 3))
        self.assertEqual(row["n_terminal_delta_aa"], -2)
        self.assertEqual(row["n_terminal_residues_a"], "MV")
        self.assertEqual(row["alignment_identity_percent"], 40)
        reverse = self.compare_row(b, a)
        self.assertEqual(reverse["protein_consequence"], row["protein_consequence"])
        self.assertEqual(reverse["n_terminal_delta_aa"], 2)
        self.assertEqual(reverse["alignment_a"], row["alignment_b"])

    def test_minus_strand_codon_mapping_preserves_biological_order(self):
        from Bio.Seq import Seq
        self.fasta.write_text(">chr\n" + str(Seq("ATGGTTGTGATGCAATAA").reverse_complement()) + "\n", encoding="utf-8")
        a = self.dataset("A", [self.cds("a", strand="-")])
        b = self.dataset("B", [self.cds("b", end=12, strand="-")])
        row = self.compare_row(a, b)
        self.assertEqual((row["alignment_a"], row["alignment_b"]), ("MVVMQ", "--MMQ"))
        self.assertEqual((row["start_delta_bp"], row["n_terminal_delta_aa"]), (6, -2))

    def test_terminal_stop_reassignment_is_not_called_a_boundary_change(self):
        self.fasta.write_text(">chr\nATGGAATGA\n", encoding="utf-8")
        a = self.dataset("A", [self.cds("a", end=9, attributes=";transl_table=4")])
        b = self.dataset("B", [self.cds("b", end=9, attributes=";transl_table=11")])
        row = self.compare_row(a, b)
        self.assertEqual(row["relationship"], "exact")
        self.assertEqual(row["stop_delta_bp"], 0)
        self.assertEqual(row["protein_consequence"], "translation_difference")
        self.assertEqual((row["alignment_a"], row["alignment_b"]), ("MEW", "ME-"))

    def test_translation_context_change_is_not_called_a_dna_mutation(self):
        self.fasta.write_text(">chr\nATGTGACAATAA\n", encoding="utf-8")
        a = self.dataset("A", [self.cds("a", end=12, attributes=";transl_table=4")])
        b = self.dataset("B", [self.cds("b", end=12, attributes=";transl_table=25")])
        row = self.compare_row(a, b)
        self.assertEqual(row["protein_consequence"], "translation_difference")
        self.assertEqual((row["internal_mismatch_count"], row["initiator_mismatch_count"]), (1, 0))
        self.assertEqual((row["alignment_a"], row["alignment_b"]), ("MWQ", "MGQ"))

    def test_phased_multipart_codons_cross_segment_boundaries_on_both_strands(self):
        for strand, genome, expected in (
            ("+", "AATGA" + "A"*15 + "AATAA", [(2, 3, 4), (5, 21, 22)]),
            ("-", "TTATT" + "A"*15 + "TCATT", [(24, 23, 22), (21, 5, 4)]),
        ):
            with self.subTest(strand=strand):
                intervals = [(1, 5), (21, 25)] if strand == "+" else [(21, 25), (1, 5)]
                features = [Feature("chr", "test", "CDS", start, end, None, strand, phase, {"ID": "joined"}) for (start, end), phase in zip(intervals, (1, 2))]
                record = next(iter_cds_sequences(features, {"chr": genome}))
                self.assertEqual(record.protein, "MK")
                self.assertEqual(cds_model(record, 25, 11, False).protein_codons(), expected)

    def test_circular_codons_wrap_without_losing_the_initiation_position(self):
        genome = "AAATAA" + "C"*9 + "ATG"
        feature = Feature("chr", "test", "CDS", 16, 24, None, "+", 0, {"ID": "wrapped"})
        record = next(iter_cds_sequences([feature], {"chr": genome}, frozenset({"chr"})))
        self.assertEqual(record.protein, "MK")
        self.assertEqual(cds_model(record, 18, 11, True).protein_codons(), [(16, 17, 18), (1, 2, 3)])

    def test_same_codon_partial_start_change_is_reported_as_initiation_context(self):
        self.fasta.write_text(">chr\nGTGAAATAA\n", encoding="utf-8")
        a = self.dataset("A", [self.cds("a", end=9)])
        b = self.dataset("B", [self.cds("b", end=9, attributes=";partial=true")])
        row = self.compare_row(a, b)
        self.assertEqual(row["relationship"], "exact")
        self.assertEqual(row["protein_consequence"], "initiator_recoding")
        self.assertEqual(row["initiator_mismatch_count"], 1)
        self.assertEqual(row["gap_aa"], 0)

    def test_parent_attributes_and_translations_are_mode_invariant(self):
        self.fasta.write_text(">chr\nGTGAAATAA\n", encoding="utf-8")
        rows = ["chr\ttest\tgene\t1\t9\t.\t+\t.\tID=parent;partial=true;gene=abc", self.cds("cds", end=9, attributes=";Parent=parent")]
        a = self.dataset("A", rows)
        ordinary = _annotation_profile(a)
        models = []
        gene = _annotation_profile(a, models)
        self.assertEqual(ordinary["gene_name_percent"], gene["gene_name_percent"])
        self.assertEqual(ordinary["gene_name_percent"], 100)
        output = self.root / "analysis"
        saved = run_analysis(a.fasta, a.gff, output, "tsv")
        self.assertEqual(saved["cds_with_gene_count"], 1)
        exported = (output / "sequences/cds_protein.fasta").read_text(encoding="utf-8")
        self.assertIn("\nVK\n", exported)
        self.assertEqual(models[0].protein, "VK")
        self.assertEqual(saved["calculation_policy"], ordinary["calculation_policy"])

    def test_cog_identifiers_do_not_require_letter_categories(self):
        a = self.dataset("A", [self.cds("a", attributes=";Dbxref=COG:COG0001")])
        b = self.dataset("B", [self.cds("b", attributes=";Dbxref=COG:COG0001")])
        profiles = [_annotation_profile(item) for item in (a, b)]
        self.assertFalse(profiles[0]["cog_data_available"])
        pair = _pairwise_rows(profiles)[0]
        self.assertEqual(pair["cog_id_jaccard"], 1)
        self.assertIsNone(pair["cog_profile_distance"])

    def test_empty_gene_value_does_not_count_as_a_name(self):
        a = self.dataset("A", [self.cds("a", attributes=";gene=")])
        self.assertEqual(_annotation_profile(a)["gene_name_percent"], 0)

    def test_zero_cds_denominators_are_unavailable(self):
        profile = _annotation_profile(self.dataset("A", []))
        for key in ("median_cds_length", "hypothetical_percent", "gene_name_percent", "recognized_start_percent"):
            self.assertIsNone(profile[key], key)

    def test_mixed_legacy_and_current_policies_are_disclosed_in_cohorts(self):
        a = self.dataset("A", [self.cds("a")])
        current = self.root / "current"
        run_analysis(a.fasta, a.gff, current, "tsv")
        legacy = self.root / "legacy"
        legacy.mkdir()
        summary = json.loads((current / "summary.json").read_text(encoding="utf-8"))
        summary.pop("calculation_policy")
        (legacy / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
        cohort = build_cohort([current, legacy])
        self.assertTrue(cohort["mixed_calculation_policies"])
        write_cohort(self.root / "cohort", cohort)
        self.assertIn("different calculation policies", (self.root / "cohort/cohort.html").read_text(encoding="utf-8"))

    def test_all_ambiguous_gc_round_trips_through_every_report(self):
        self.fasta.write_text(">chr\nNNNNNNNNN\n", encoding="utf-8")
        a, b = self.dataset("A", []), self.dataset("B", [])
        output = self.root / "unknown"
        saved = run_analysis(a.fasta, a.gff, output, "tsv")
        self.assertIsNone(saved["quality_control"]["genome_gc_percent"])
        self.assertIn("Not available", (output / "report.html").read_text(encoding="utf-8"))
        cohort = build_cohort([output])
        self.assertIsNone(cohort["samples"][0]["genome_gc_percent"])
        write_cohort(self.root / "cohort", cohort)
        compared = run_comparison([a, b], self.root / "comparison", gene_level=True)
        self.assertIsNone(compared["datasets"][0]["gc_percent"])
        self.assertIsNone(compared["pairwise_comparisons"][0]["gc_percent_delta"])
        command = [sys.executable, "-m", "annostat", "analyze", "-f", str(a.fasta), "-g", str(a.gff), "-o", str(self.root / "cli")]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("not available", result.stdout.lower())
        self.fasta.write_text(">chr\nATATATATA\n", encoding="utf-8")
        self.assertEqual(_annotation_profile(a)["gc_percent"], 0)


if __name__ == "__main__":
    unittest.main()
