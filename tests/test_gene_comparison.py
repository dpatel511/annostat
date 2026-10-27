"""Gene-level comparison of two annotations on unchanged genome sequences."""

from __future__ import annotations

import csv
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from annostat.comparison import GenomeInput, run_comparison
from annostat.gene_comparison import CdsModel, _candidate_links
from annostat.models import Feature


class GeneComparisonTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.fasta = self.root / "genome.fna"
        self.fasta.write_text(">chr\n" + "ATG" * 100 + "\n", encoding="utf-8")

    def dataset(self, label, rows, fasta=None):
        gff = self.root / f"{label}.gff3"
        gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return GenomeInput(label, fasta or self.fasta, gff)

    def cds(self, identifier, start, end, strand="+", phase=0, attributes=""):
        return f"chr\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID={identifier}{attributes}"

    def compare(self, left, right):
        return run_comparison([left, right], self.root / "results", gene_level=True)

    def test_exact_models_ignore_ids_and_attribute_order(self):
        a = self.dataset("A", [self.cds("old", 1, 30, attributes=";gene=abc;Dbxref=COG:COG0001,GO:1")])
        b = self.dataset("B", [self.cds("new", 1, 30, attributes=";Dbxref=GO:1,COG:COG0001;gene=abc")])
        summary = self.compare(a, b)
        result = summary["gene_level_comparison"]
        row = result["relationships"][0]
        self.assertEqual((row["status"], row["relationship"]), ("matched", "exact"))
        self.assertEqual(row["changed_attributes"], "")
        self.assertTrue(row["protein_identical"])
        self.assertEqual(result["exact_agreement"], 1.0)
        self.assertEqual(result["coordinate_system"], "1-based inclusive")

    def test_biological_boundary_deltas_on_both_strands(self):
        for strand, start, end in (("+", 4, 30), ("-", 1, 27)):
            with self.subTest(strand=strand):
                a = self.dataset("A", [self.cds("a", 1, 30, strand)])
                b = self.dataset("B", [self.cds("b", start, end, strand)])
                row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
                self.assertEqual(row["relationship"], "start_shift")
                self.assertEqual(row["start_delta_bp"], 3)
                self.assertEqual(row["stop_delta_bp"], 0)
                self.assertEqual(row["cds_length_b"] - row["cds_length_a"], -3)
                reverse = self.compare(b, a)["gene_level_comparison"]["relationships"][0]
                self.assertEqual(reverse["start_delta_bp"], -3)

    def test_frame_disagreement_is_not_called_a_start_shift(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 2, 30)])
        row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
        self.assertEqual(row["relationship"], "frame_disagreement")
        self.assertEqual(row["status"], "candidate")

    def test_exact_anchors_do_not_turn_real_overlaps_into_ambiguity(self):
        a = self.dataset("A", [self.cds("a1", 1, 30), self.cds("a2", 25, 54)])
        b = self.dataset("B", [self.cds("b2", 25, 54), self.cds("b1", 1, 30)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["status_counts"], {"matched": 2})
        self.assertEqual(result["exact_agreement"], 1.0)

    def test_many_to_one_overlaps_remain_ambiguous_and_unmatched_are_retained(self):
        a = self.dataset("A", [self.cds("a1", 1, 15), self.cds("a2", 16, 30), self.cds("only_a", 61, 90)])
        b = self.dataset("B", [self.cds("joined", 1, 30), self.cds("only_b", 121, 150)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["status_counts"], {"ambiguous": 2, "unmatched": 2})
        self.assertEqual({row["id_a"] for row in result["relationships"] if row["id_a"]}, {"a1", "a2", "only_a"})
        self.assertEqual({row["id_b"] for row in result["relationships"] if row["id_b"]}, {"joined", "only_b"})

    def test_duplicate_models_have_exact_agreement_without_inventing_id_matches(self):
        a = self.dataset("A", [self.cds("a1", 1, 30), self.cds("a2", 1, 30)])
        result = self.compare(a, a.__class__("B", a.fasta, a.gff))["gene_level_comparison"]
        self.assertEqual(result["exact_agreement"], 1.0)
        self.assertEqual(result["status_counts"], {"ambiguous": 4})

    def test_ambiguous_exact_models_are_not_overruled_by_a_boundary_alternative(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b1", 1, 30), self.cds("b2", 1, 30), self.cds("alternative", 4, 30)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["status_counts"], {"ambiguous": 3})
        self.assertEqual(result["exact_model_count"], 1)

    def test_multipart_cds_are_counted_once_and_gaps_are_not_overlap(self):
        a = self.dataset("A", [self.cds("joined", 1, 9), self.cds("joined", 31, 39)])
        b = self.dataset("B", [self.cds("gap", 16, 24)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["cds_counts"], {"A": 1, "B": 1})
        self.assertEqual(result["status_counts"], {"unmatched": 2})

    def test_annotation_changes_include_parent_pseudo_and_translation_table(self):
        parent = "chr\ttest\tgene\t1\t30\t.\t+\t.\tID=parent;pseudo=true"
        a = self.dataset("A", [parent, self.cds("a", 1, 30, attributes=";Parent=parent;product=hypothetical protein;transl_table=11")])
        b = self.dataset("B", [self.cds("b", 1, 30, attributes=";product=named protein;transl_table=4")])
        row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
        self.assertEqual(row["relationship"], "exact")
        self.assertTrue({"product", "pseudo", "genetic_code"} <= set(row["changed_attributes"].split(";")))

    def test_fasta_identity_checks_bases_and_sequence_ids_not_formatting(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        other = self.root / "other.fna"
        other.write_text(">chr description\n" + "atg" * 50 + "\n" + "atg" * 50 + "\n", encoding="utf-8")
        b = self.dataset("B", [self.cds("b", 1, 30)], other)
        self.assertEqual(self.compare(a, b)["gene_level_comparison"]["exact_agreement"], 1.0)
        other.write_text(">chr\n" + "ATG" * 99 + "ATA\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "identical genome sequences"):
            self.compare(a, b)

    def test_cli_reports_gene_results_and_default_mode_has_no_gene_output(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 4, 30, attributes=";product=changed")])
        command = [sys.executable, "-m", "annostat", "compare", "--genome", "A", str(a.fasta), str(a.gff), "--genome", "B", str(b.fasta), str(b.gff)]
        default = self.root / "default"
        run_comparison([a, b], default)
        self.assertFalse((default / "tables" / "gene_comparisons.tsv").exists())
        output = self.root / "cli"
        completed = subprocess.run(command + ["--gene-level", "-o", str(output)], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("Relationship rows", completed.stdout)
        self.assertIn("Accepted CDS pairs", completed.stdout)
        saved = json.loads((output / "comparison.json").read_text(encoding="utf-8"))
        with (output / "tables" / "gene_comparisons.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), len(saved["gene_level_comparison"]["relationships"]))
        report = (output / "comparison.html").read_text(encoding="utf-8")
        self.assertIn("CDS annotation comparison", report)
        self.assertIn("2 annotations", report)
        self.assertNotIn("2 assemblies", report)
        self.assertIn("gene_comparisons.tsv", report)
        self.assertIn("1-based inclusive", report)

    def test_gene_level_requires_exactly_two_annotations(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        with self.assertRaisesRegex(ValueError, "exactly two"):
            run_comparison([a, a.__class__("B", a.fasta, a.gff), a.__class__("C", a.fasta, a.gff)], self.root / "three", gene_level=True)

    def test_stop_shifts_follow_strand_direction(self):
        for strand, start, end in (("+", 1, 33), ("-", 1, 33)):
            with self.subTest(strand=strand):
                a_start, a_end = (1, 30) if strand == "+" else (4, 33)
                a = self.dataset("A", [self.cds("a", a_start, a_end, strand)])
                b = self.dataset("B", [self.cds("b", start, end, strand)])
                row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
                self.assertEqual(row["relationship"], "stop_shift")
                self.assertEqual((row["start_delta_bp"], row["stop_delta_bp"]), (0, 3))

    def test_same_coordinates_with_different_phase_are_not_exact(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 1, 30, phase=1)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["exact_agreement"], 0)
        self.assertEqual(result["relationships"][0]["relationship"], "frame_disagreement")

    def test_minus_multipart_models_are_order_independent(self):
        a = self.dataset("A", [self.cds("a", 1, 9, "-"), self.cds("a", 31, 39, "-")])
        b = self.dataset("B", [self.cds("b", 31, 39, "-"), self.cds("b", 1, 9, "-")])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["exact_agreement"], 1)
        self.assertEqual(result["relationships"][0]["cds_length_a"], 18)
        self.assertEqual(result["relationships"][0]["segments_a"], "31-39;1-9")

    def test_circular_boundaries_use_biological_order(self):
        region = "chr\ttest\tregion\t1\t300\t.\t+\t.\tID=chr;Is_circular=true"
        for strand, start, end in (("+", 294, 315), ("-", 291, 312)):
            with self.subTest(strand=strand):
                a = self.dataset("A", [region, self.cds("a", 291, 315, strand)])
                b = self.dataset("B", [region, self.cds("b", start, end, strand)])
                row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
                self.assertEqual(row["relationship"], "start_shift")
                self.assertEqual((row["start_delta_bp"], row["stop_delta_bp"]), (3, 0))
                self.assertEqual(row["overlap_bp"], 22)

    def test_boundary_alternatives_can_cross_the_circular_origin(self):
        region = "chr\ttest\tregion\t1\t300\t.\t+\t.\tID=chr;Is_circular=true"
        a = self.dataset("A", [region, self.cds("a", 298, 315)])
        b = self.dataset("B", [region, self.cds("b", 1, 15)])
        row = self.compare(a, b)["gene_level_comparison"]["relationships"][0]
        self.assertEqual(row["relationship"], "start_shift")
        self.assertEqual((row["start_delta_bp"], row["stop_delta_bp"]), (3, 0))
        swapped = self.compare(b, a)["gene_level_comparison"]["relationships"][0]
        self.assertEqual(swapped["start_delta_bp"], -3)

    def test_full_circle_models_preserve_their_start_position(self):
        region = "chr\ttest\tregion\t1\t300\t.\t+\t.\tID=chr;Is_circular=true"
        a = self.dataset("A", [region, self.cds("a", 1, 300)])
        b = self.dataset("B", [region, self.cds("b", 4, 303)])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertEqual(result["exact_agreement"], 0)
        self.assertEqual(result["relationships"][0]["overlap_bp"], 300)

    def test_empty_cds_sets_write_headers_and_have_no_agreement_denominator(self):
        a = self.dataset("A", [])
        b = self.dataset("B", [])
        result = self.compare(a, b)["gene_level_comparison"]
        self.assertIsNone(result["exact_agreement"])
        self.assertEqual(result["relationships"], [])
        with (self.root / "results" / "tables" / "gene_comparisons.tsv").open(encoding="utf-8", newline="") as handle:
            table = csv.DictReader(handle, delimiter="\t")
            self.assertIn("id_a", table.fieldnames)
            self.assertEqual(list(table), [])

    def test_unknown_cds_strands_are_rejected_and_opposite_strands_do_not_match(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 1, 30, "-")])
        self.assertEqual(self.compare(a, b)["gene_level_comparison"]["status_counts"], {"unmatched": 2})
        b = self.dataset("B", [self.cds("b", 1, 30, ".")])
        with self.assertRaisesRegex(ValueError, "explicit.*strands"):
            self.compare(a, b)

    def test_default_rerun_removes_stale_gene_table_and_report_section(self):
        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 1, 30)])
        self.compare(a, b)
        output = self.root / "results"
        summary = run_comparison([a, b], output)
        self.assertNotIn("gene_level_comparison", summary)
        self.assertFalse((output / "tables" / "gene_comparisons.tsv").exists())
        self.assertNotIn("CDS annotation comparison", (output / "comparison.html").read_text(encoding="utf-8"))

    def test_interval_sweep_matches_independent_base_set_oracle(self):
        rng = random.Random(42)
        groups = []
        for side in range(2):
            models = []
            for index in range(50):
                start = rng.randrange(1, 200)
                end = start + rng.randrange(1, 40)
                strand = rng.choice(("+", "-"))
                seqid = rng.choice(("chr", "plasmid"))
                blocks = [(start, end, 0)]
                if index % 3 == 0:
                    blocks.append((end + 10, end + 20, 0))
                feature = Feature(seqid, "test", "CDS", start, end, None, strand, 0, {"ID": f"{side}-{index}"})
                models.append(CdsModel(feature, tuple(blocks), end - start + 1, 0, "", 11, False, 300))
            groups.append(models)
        left, right = _candidate_links(*groups)

        def bases(model):
            return {position for start, end, _ in model.blocks for position in range(start, end + 1)}

        expected = {
            (i, j)
            for i, a in enumerate(groups[0])
            for j, b in enumerate(groups[1])
            if (a.feature.seqid, a.feature.strand) == (b.feature.seqid, b.feature.strand)
            and bases(a) & bases(b)
        }
        self.assertEqual({(i, j) for i, neighbours in enumerate(left) for j in neighbours}, expected)
        self.assertEqual({(i, j) for j, neighbours in enumerate(right) for i in neighbours}, expected)

    def test_row_reordering_preserves_relationships_and_swapping_preserves_counts(self):
        a_rows = [self.cds("a1", 1, 30), self.cds("a2", 61, 90)]
        b_rows = [self.cds("b1", 4, 30), self.cds("b2", 61, 90)]
        a, b = self.dataset("A", a_rows), self.dataset("B", b_rows)
        first = self.compare(a, b)["gene_level_comparison"]
        reordered = self.compare(self.dataset("A", a_rows[::-1]), self.dataset("B", b_rows[::-1]))["gene_level_comparison"]
        self.assertEqual(first, reordered)
        swapped = self.compare(b, a)["gene_level_comparison"]
        self.assertEqual(first["exact_agreement"], swapped["exact_agreement"])
        self.assertEqual(first["status_counts"], swapped["status_counts"])

    def test_gene_comparison_does_not_reparse_or_retranslate_inputs(self):
        from unittest.mock import patch
        from annostat.parsers import parse_fasta, parse_gff
        from annostat.sequences import iter_cds_sequences

        a = self.dataset("A", [self.cds("a", 1, 30)])
        b = self.dataset("B", [self.cds("b", 4, 30)])
        with patch("annostat.validation.parse_fasta", wraps=parse_fasta) as fasta_parser, patch("annostat.validation.parse_gff", wraps=parse_gff) as gff_parser, patch("annostat.comparison.iter_cds_sequences", wraps=iter_cds_sequences) as cds_iterator:
            self.compare(a, b)
        self.assertEqual(fasta_parser.call_count, 2)
        self.assertEqual(gff_parser.call_count, 2)
        self.assertEqual(cds_iterator.call_count, 2)


if __name__ == "__main__":
    unittest.main()
