"""Unapplied source exceptions must not become definitive protein evidence."""
import unittest
from dataclasses import replace
from annostat.models import Feature
from annostat.sequences import iter_cds_sequences
from annostat.qc import sequence_quality_findings
from annostat.gene_comparison import cds_model, compare_cds


class UnappliedExceptionTests(unittest.TestCase):
    def record(self, annotation, *, strand="+", phase=0):
        from Bio.Seq import Seq
        dna = "A" * phase + "ATGAAAGGGTAA"
        genome = str(Seq(dna).reverse_complement()) if strand == "-" else dna
        feature = Feature("chr", "test", "CDS", 1, len(genome), None, strand, phase,
                          {"ID": "cds", "transl_except": annotation})
        return next(iter(iter_cds_sequences([feature], {"chr": genome})))

    def test_unsupported_unknown_and_unmapped_exceptions_are_visible(self):
        for annotation in ("(pos:join(4..5,6),aa:Sec)", "(pos:4..6,aa:MadeUp)",
                           "(pos:2..4,aa:Sec)", "(pos:100..102,aa:Sec)",
                           "(pos:6..4,aa:Sec)", "broken"):
            with self.subTest(annotation=annotation):
                record = self.record(annotation)
                self.assertEqual(record.unapplied_translation_exceptions, (annotation,))
                findings = sequence_quality_findings(record)
                issue = next(x for x in findings if x.issue_type == "unapplied_translation_exception")
                self.assertEqual(issue.severity, "warning")
                self.assertIn(annotation, issue.details)

    def test_valid_and_invalid_entries_are_accounted_separately(self):
        annotation = "(pos:4..6,aa:Sec),(pos:join(7..8,9),aa:Leu)"
        record = self.record(annotation)
        self.assertEqual(record.protein, "MUG")
        self.assertEqual(record.translation_exception_indices, (1,))
        self.assertEqual(record.unapplied_translation_exceptions, ("(pos:join(7..8,9),aa:Leu)",))

    def test_supported_entries_and_native_terminal_term_are_not_flagged(self):
        for strand, annotation in (("+", "(pos:4..6,aa:Sec),(pos:10..12,aa:TERM)"),
                                   ("-", "(pos:complement(7..9),aa:Sec),(pos:complement(1..3),aa:TERM)")):
            record = self.record(annotation, strand=strand)
            self.assertEqual(record.protein, "MUG")
            self.assertEqual(record.unapplied_translation_exceptions, ())
            self.assertNotIn("unapplied_translation_exception", [x.issue_type for x in sequence_quality_findings(record)])

    def test_phase_trimmed_exception_is_unapplied(self):
        record = self.record("(pos:1..3,aa:Sec)", phase=1)
        self.assertTrue(record.unapplied_translation_exceptions)

    def test_unapplied_exception_forces_limited_evidence_in_both_directions(self):
        record = self.record("(pos:join(4..5,6),aa:Leu)")
        other = replace(record, feature=replace(record.feature, attributes={"ID": "other"}),
                        unapplied_translation_exceptions=())
        a, b = cds_model(record, 12, 11, False), cds_model(other, 12, 11, False)
        for left, right in ((a, b), (b, a)):
            row = compare_cds([left], [right], ("A", "B"))["relationships"][0]
            self.assertEqual(row["status"], "matched")
            self.assertTrue(row["protein_identical"])
            self.assertEqual(row["protein_consequence"], "limited")
            self.assertIsNone(row["alignment_identity_percent"])
            self.assertIsNone(row["substitution_count"])
            self.assertIn("Unapplied", row["consequence_note"])
            self.assertIn("join(4..5,6)", row["unapplied_transl_except_a"] + row["unapplied_transl_except_b"])

    def test_multipart_repeated_exception_is_not_duplicated(self):
        annotation = "(pos:join(4..5,10),aa:Leu)"
        parts = [Feature("chr", "test", "CDS", 1, 5, None, "+", 0, {"ID": "x", "transl_except": annotation}),
                 Feature("chr", "test", "CDS", 10, 16, None, "+", 1, {"ID": "x", "transl_except": annotation})]
        record = next(iter(iter_cds_sequences(parts, {"chr": "ATGAACCCCAAAGTAA"})))
        self.assertEqual(record.unapplied_translation_exceptions, (annotation,))

    def test_circular_unwrapped_exception_still_maps_to_its_codon(self):
        # 1-based inclusive: ATG at 7,8,9; AAA at 10,11,12; TAA at 1,2,3.
        feature = Feature("chr", "test", "CDS", 7, 15, None, "+", 0,
                          {"ID": "wrap", "transl_except": "(pos:10..12,aa:Sec)"})
        record = next(iter(iter_cds_sequences([feature], {"chr": "TAACCCATGAAA"}, frozenset({"chr"}))))
        self.assertEqual(record.protein, "MU")
        self.assertEqual(record.unapplied_translation_exceptions, ())

    def test_public_check_rejects_a_changed_snapshot_before_running_annostat(self):
        import subprocess
        import sys
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "bakta.fna").write_text(">changed\nATG\n")
            script = Path(__file__).resolve().parents[1] / "examples/check_public_comparison.py"
            result = subprocess.run([sys.executable, str(script), "--inputs", str(root),
                                     "--output", str(root / "result")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("Frozen source mismatch", result.stderr)
            self.assertFalse((root / "result").exists())
