"""Description differences must not be mistaken for protein or function changes."""
import csv
import json
import re
import tempfile
import unittest
from pathlib import Path
from urllib.parse import quote

from annostat.comparison import GenomeInput, run_comparison
from annostat.gene_comparison import product_status
from annostat.gene_report import _description_percent


class DescriptionSummaryTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.fasta = self.root / "genome.fna"
        self.fasta.write_text(">chr\n" + "ATG" * 200 + "\n", encoding="utf-8")

    def dataset(self, label, models):
        gff = self.root / (label + ".gff3")
        rows = []
        for identifier, start, end, product, gene in models:
            attrs = "ID=" + identifier
            if product is not None:
                attrs += ";product=" + quote(product)
            if gene is not None:
                attrs += ";gene=" + quote(gene)
            rows.append(f"chr\ttest\tCDS\t{start}\t{end}\t.\t+\t0\t{attrs}")
        gff.write_text("##gff-version 3\n" + "\n".join(rows) + "\n", encoding="utf-8")
        return GenomeInput(label, self.fasta, gff)

    def matrix_inputs(self):
        products = [None, "HYPOTHETICAL protein", "DUF1234 domain-containing protein"]
        choices = [(a, b) for a in range(3) for b in range(3)] + [(2, 0)]
        left, right = [], []
        for i, (a, b) in enumerate(choices):
            start = 1 + i * 36  # GFF3 coordinates are 1-based inclusive.
            left.append((f"a{i}", start, start + 29, products[a], "abc" if i % 2 == 0 else None))
            right.append((f"b{i}", start, start + 29, products[b], "xyz" if i < 9 and i % 2 else None))
        return [self.dataset("SourceA", left), self.dataset("SourceB", right)]

    def test_small_nonzero_percentages_are_not_displayed_as_zero(self):
        self.assertEqual(_description_percent(100 * 2 / 4213), "<0.1%")
        self.assertEqual(_description_percent(0), "0.0%")
        self.assertEqual(_description_percent(None), "Not available")
        self.assertEqual(_description_percent(30.4), "30.4%")

    def test_product_status_is_a_label_convention_not_function_validation(self):
        for value, expected in [
            (None, "missing"), ("   ", "missing"),
            ("conserved HYPOTHETICAL protein", "explicit_hypothetical"),
            ("DUF1234 protein", "other_product_text"),
            ("unknown protein", "other_product_text"),
        ]:
            with self.subTest(product=value):
                self.assertEqual(product_status(value), expected)

    def test_directional_matrix_reconciles_exports_and_reverses_with_source_order(self):
        inputs = self.matrix_inputs()
        out = self.root / "forward"
        summary = run_comparison(inputs, out, gene_level=True)
        gene = summary["gene_level_comparison"]
        description = gene["description_summary"]
        self.assertEqual(gene["schema_version"], "1.7")
        self.assertEqual(description["identical_protein_pairs"], 10)
        self.assertEqual(description["product_text_difference_pairs"], 7)
        self.assertEqual(description["unresolved_cds_counts"], {"SourceA": 0, "SourceB": 0})
        self.assertEqual(description["gene_symbol_counts"], {
            "SourceA": {"provided": 5, "missing": 5},
            "SourceB": {"provided": 4, "missing": 6},
        })
        matrix = description["product_status_matrix"]
        self.assertEqual(len(matrix), 9)
        for cell in matrix:
            expected = 2 if (cell["product_status_a"], cell["product_status_b"]) == ("other_product_text", "missing") else 1
            self.assertEqual(cell["pair_count"], expected)
            self.assertEqual(cell["percent_of_identical_pairs"], 10 * expected)
        self.assertEqual(sum(cell["pair_count"] for cell in matrix), 10)
        reverse = run_comparison(inputs[::-1], self.root / "reverse", gene_level=True)["gene_level_comparison"]["description_summary"]
        transposed = {(cell["product_status_b"], cell["product_status_a"]): cell["pair_count"] for cell in reverse["product_status_matrix"]}
        self.assertEqual({(cell["product_status_a"], cell["product_status_b"]): cell["pair_count"] for cell in matrix}, transposed)
        self.assertEqual(reverse["gene_symbol_counts"], description["gene_symbol_counts"])
        with (out / "tables/gene_product_status.tsv").open(encoding="utf-8", newline="") as handle:
            table = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(table), 9)
        for saved, cell in zip(table, matrix):
            self.assertEqual(saved["dataset_a"], "SourceA")
            self.assertEqual(saved["dataset_b"], "SourceB")
            self.assertEqual(int(saved["identical_protein_pairs"]), 10)
            self.assertEqual(int(saved["pair_count"]), cell["pair_count"])
            self.assertEqual(float(saved["percent_of_identical_pairs"]), cell["percent_of_identical_pairs"])
        with (out / "tables/gene_comparisons.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["product_status_a"], "missing")
        self.assertEqual(rows[0]["product_status_b"], "missing")
        html = (out / "comparison.html").read_text(encoding="utf-8")
        payload = json.loads(re.search(r'<script id="gene-data" type="application/json">(.*?)</script>', html, re.S).group(1))
        self.assertEqual(payload["description_summary"], description)
        self.assertIn("of 10 identical-protein pairs differ in product text", html)
        self.assertIn('data-product-a="other_product_text" data-product-b="missing"', html)
        self.assertIn("Other product text does not establish a known or experimentally characterized function", html)
        self.assertIn("<b>7</b> of 10 identical-protein pairs differ in product text", html)
        text = (out / "comparison.txt").read_text(encoding="utf-8")
        self.assertIn("Product-label comparison: 10 accepted pairs with reliable identical predicted proteins", text)
        self.assertIn("A Other product text / B Missing product: 2 / 10 (20.0%)", text)
        self.assertIn("SourceA: gene symbol provided 5 / 10; missing 5", text)
        self.assertIn("SourceB: gene symbol provided 4 / 10; missing 6", text)
        run_comparison(inputs, out)
        self.assertFalse((out / "tables/gene_product_status.tsv").exists())

    def test_limited_equal_strings_and_many_to_many_rows_stay_out_of_matrix(self):
        self.fasta.write_text(">chr\nNNN" + "ATG" * 199 + "\n", encoding="utf-8")
        left = [("limited_a", 1, 30), ("changed_a", 61, 90), ("a1", 121, 150), ("a2", 151, 180)]
        right = [("limited_b", 1, 30), ("changed_b", 64, 90), ("b1", 136, 165), ("b2", 166, 195)]
        inputs = [self.dataset(label, [(identifier, start, end, "same description", None) for identifier, start, end in models])
                  for label, models in (("A", left), ("B", right))]
        out = self.root / "limited"
        gene = run_comparison(inputs, out, gene_level=True)["gene_level_comparison"]
        limited = next(row for row in gene["relationships"] if row["id_a"] == "limited_a")
        self.assertTrue(limited["protein_identical"])
        self.assertEqual(limited["protein_consequence"], "limited")
        self.assertEqual(gene["status_counts"]["ambiguous"], 3)
        description = gene["description_summary"]
        self.assertEqual(description["identical_protein_pairs"], 0)
        self.assertEqual(description["interpretable_protein_difference_pairs"], 1)
        self.assertEqual(description["limited_protein_pairs"], 1)
        self.assertEqual(description["unresolved_cds_counts"], {"A": 2, "B": 2})
        self.assertTrue(all(cell["pair_count"] == 0 and cell["percent_of_identical_pairs"] is None
                            for cell in description["product_status_matrix"]))
        self.assertIn("Unresolved CDS models (distinct): A: 2; B: 2", (out / "comparison.txt").read_text(encoding="utf-8"))
        with (out / "tables/gene_product_status.tsv").open(encoding="utf-8", newline="") as handle:
            self.assertTrue(all(row["percent_of_identical_pairs"] == "" for row in csv.DictReader(handle, delimiter="\t")))
        html = (out / "comparison.html").read_text(encoding="utf-8")
        self.assertIn("No accepted pairs with reliable identical predicted proteins", html)

    def test_empty_annotations_and_blank_products_preserve_missingness(self):
        empty = [self.dataset(label, []) for label in ("A", "B")]
        description = run_comparison(empty, self.root / "empty", gene_level=True)["gene_level_comparison"]["description_summary"]
        self.assertEqual(description["accepted_pairs"], 0)
        self.assertEqual(description["gene_symbol_counts"], {"A": {"provided": 0, "missing": 0}, "B": {"provided": 0, "missing": 0}})
        blanks = [self.dataset("A", [("a", 1, 30, "  ", " "), ("only_a", 181, 210, "DUF1234 protein", None)]),
                  self.dataset("B", [("b", 1, 30, None, None)])]
        gene = run_comparison(blanks, self.root / "blank", gene_level=True)["gene_level_comparison"]
        self.assertEqual(gene["description_summary"]["product_status_matrix"][0]["pair_count"], 1)
        self.assertEqual(gene["description_summary"]["gene_symbol_counts"]["A"], {"provided": 0, "missing": 1})
        self.assertEqual(gene["description_summary"]["product_text_difference_pairs"], 0)
        unmatched = next(row for row in gene["relationships"] if row["id_a"] == "only_a")
        self.assertEqual(unmatched["product_status_a"], "other_product_text")
        self.assertIsNone(unmatched["product_status_b"])


if __name__ == "__main__":
    unittest.main()
