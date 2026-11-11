"""Source-context facts are separate from accepted CDS correspondence (1-based inclusive)."""
import csv
import json
import re
import tempfile
import unittest
from pathlib import Path

from annostat.comparison import GenomeInput, run_comparison
from annostat.gene_comparison import CdsModel, compare_cds
from annostat.models import Feature
from annostat.locus_context import collect_source_context


def feature(kind, start, end, identifier, strand="+", seqid="chr", **attributes):
    return Feature(seqid, "test", kind, start, end, None, strand,
                   0 if kind == "CDS" else None, dict(ID=identifier, **attributes))


def model(item, blocks=None, circular=False, length=200):
    blocks = blocks or ((item.start, item.end, 0),)
    return CdsModel(item, blocks, sum(end-start+1 for start, end, _ in blocks),
                    1, "", 11, circular, length, "M")


class SourceContextTests(unittest.TestCase):
    def context(self, a, b, features_a=(), features_b=()):
        gene = compare_cds(a, b, ("A", "B"))
        before = json.dumps(gene, sort_keys=True)
        result = collect_source_context(gene, (a, b), (features_a, features_b), ("A", "B"))
        self.assertEqual(json.dumps(gene, sort_keys=True), before)
        return result

    def test_overlap_facts_keep_source_ids_strands_and_explicit_denominators(self):
        a = model(feature("CDS", 10, 40, "a"))
        pseudo = feature("gene", 20, 49, "p", pseudogene="unknown")
        rna = feature("tRNA", 30, 50, "r", strand="-")
        opposite = feature("CDS", 25, 54, "o", strand="-")
        result = self.context([a], [model(opposite)], [a.feature], [pseudo, rna, opposite])
        facts = {r["observation"]: r for r in result["observations"] if r["dataset"] == "A"}
        self.assertEqual(set(facts), {"gene_only_pseudogene_overlap", "rna_overlap", "opposite_strand_cds_overlap"})
        for kind, identifier, bp, total in (("gene_only_pseudogene_overlap", "p", 21, 30),
                                           ("rna_overlap", "r", 11, 21),
                                           ("opposite_strand_cds_overlap", "o", 16, 30)):
            row = facts[kind]
            self.assertEqual((row["cds_id"], row["context_dataset"], row["feature_id"]), ("a", "B", identifier))
            self.assertEqual((row["overlap_bp"], row["cds_covered_bp"], row["feature_covered_bp"]), (bp, 31, total))
            self.assertAlmostEqual(row["cds_coverage_percent"], 100*bp/31)
            self.assertAlmostEqual(row["feature_coverage_percent"], 100*bp/total)
        self.assertEqual(facts["rna_overlap"]["feature_strand"], "-")

    def test_gene_only_excludes_cds_descendants_but_not_rna_and_false_is_not_true(self):
        a = model(feature("CDS", 10, 40, "a"))
        rows = [feature("gene", 10, 40, "parent", pseudo="true"),
                feature("mRNA", 10, 40, "transcript", Parent="parent"),
                feature("CDS", 10, 40, "child", Parent="transcript"),
                feature("gene", 10, 40, "not_pseudo", pseudo="false"),
                feature("pseudogene", 10, 40, "gene_only"),
                feature("ncRNA", 10, 40, "rna", Parent="gene_only"),
                feature("tRNA", 10, 40, "elsewhere", seqid="other")]
        result = self.context([a], [], [a.feature], rows)
        self.assertEqual({r["feature_id"] for r in result["observations"]}, {"gene_only", "rna"})

    def test_circular_and_repeated_blocks_use_physical_union_not_enclosing_span(self):
        a = model(feature("CDS", 95, 110, "a"), ((95, 100, 0), (1, 10, 0), (5, 10, 0)),
                  circular=True, length=100)
        rna = feature("rRNA", 98, 105, "r")
        rows = [feature("region", 1, 100, "chr", Is_circular="true"), rna,
                feature("tRNA", 40, 50, "gap")]
        result = self.context([a], [], [a.feature], rows)
        fact = next(r for r in result["observations"] if r["observation"] == "rna_overlap")
        self.assertEqual((fact["overlap_bp"], fact["cds_covered_bp"], fact["feature_covered_bp"]), (8, 16, 8))
        self.assertEqual(fact["cds_coverage_percent"], 50)
        self.assertEqual(fact["feature_coverage_percent"], 100)
        self.assertEqual(fact["feature_segments"], "1-5;98-100")
        self.assertEqual(json.loads(fact["source_attributes"])[0]["end"], 105)
        self.assertNotIn("gap", {r["feature_id"] for r in result["observations"]})

    def test_sequence_edge_does_not_invent_circularity_or_truncation(self):
        a = model(feature("CDS", 1, 9, "a"))
        result = self.context([a], [], [a.feature], [])
        row = result["observations"][0]
        self.assertEqual(row["observation"], "sequence_boundary")
        self.assertIsNone(row["overlap_bp"])
        self.assertIn("not declared", row["note"])
        self.assertIn("does not establish truncation", row["note"])
        a = model(a.feature, circular=True)
        row = self.context([a], [], [a.feature], [])["observations"][0]
        self.assertIn("circular", row["note"])

    def test_one_context_per_distinct_model_not_per_ambiguous_link(self):
        a = [model(feature("CDS", 10, 39, "a")), model(feature("CDS", 40, 69, "b"))]
        b = [model(feature("CDS", 20, 59, "x"))]
        rna = feature("tRNA", 10, 69, "r")
        result = self.context(a, b, [m.feature for m in a], [b[0].feature, rna])
        facts = [r for r in result["observations"] if r["observation"] == "rna_overlap"]
        self.assertEqual([r["cds_id"] for r in facts], ["a", "b"])
        self.assertEqual(result["model_count"], 3)
        self.assertEqual(result["models_with_observations"], 2)

    def test_accepted_pairs_and_empty_inputs_have_no_unresolved_context(self):
        a = model(feature("CDS", 10, 39, "a"))
        b = model(feature("CDS", 10, 39, "b"))
        result = self.context([a], [b], [a.feature], [b.feature, feature("rRNA", 10, 39, "r")])
        self.assertEqual(result["observations"], [])
        self.assertEqual(result["model_count"], 0)
        self.assertEqual(self.context([], [])["observations"], [])

    def test_multipart_rna_counts_once_and_preserves_rows(self):
        a = model(feature("CDS", 10, 40, "a"))
        parts = [feature("tRNA", 10, 20, "r"), feature("tRNA", 18, 25, "r")]
        result = self.context([a], [], [a.feature], parts)
        row = result["observations"][0]
        self.assertEqual((row["overlap_bp"], row["feature_covered_bp"]), (16, 16))
        self.assertEqual(len(json.loads(row["source_attributes"])), 2)

    def test_sweep_matches_independent_base_set_oracle(self):
        import random
        rng = random.Random(91)
        a = [model(feature("CDS", start, start+rng.randrange(3, 20), f"a{i}"))
             for i, start in enumerate(rng.sample(range(5, 170), 35))]
        rnas = [feature("ncRNA", start, start+rng.randrange(1, 25), f"r{i}", strand=".")
                for i, start in enumerate(rng.sample(range(2, 170), 40))]
        result = self.context(a, [], [m.feature for m in a], rnas)
        actual = {(r["cds_id"], r["feature_id"]): r["overlap_bp"] for r in result["observations"]}
        expected = {}
        for m in a:
            for r in rnas:
                overlap = len(set(range(m.feature.start, m.feature.end+1)) & set(range(r.start, r.end+1)))
                if overlap:
                    expected[m.feature.id, r.id] = overlap
        self.assertEqual(actual, expected)


class SourceContextExportTests(unittest.TestCase):
    def test_json_tsv_html_and_text_round_trip_source_facts_and_aggregate_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "genome.fna"
            fasta.write_text(">chr\n" + "ATG"*40 + "\n")
            a, b = root / "a.gff3", root / "b.gff3"
            a.write_text("##gff-version 3\nchr\ttest\tCDS\t10\t39\t.\t+\t0\tID=a;partial=true\n")
            b.write_text("##gff-version 3\nchr\ttest\tgene\t20\t49\t.\t-\t.\tID=p;pseudo=true;Note=%3C/script%3E%3Cimg%20src=x%20onerror=alert(1)%3E\n")
            inputs = [GenomeInput("A", fasta, a), GenomeInput("B", fasta, b)]
            out = root / "out"
            summary = run_comparison(inputs, out, gene_level=True)
            gene = summary["gene_level_comparison"]
            self.assertEqual(gene["schema_version"], "1.5")
            observation = gene["source_context"]["observations"][0]
            self.assertEqual((observation["overlap_bp"], observation["feature_strand"]), (20, "-"))
            self.assertEqual(gene["status_counts"], {"unmatched": 1})
            with (out / "tables/gene_context.tsv").open(newline="") as handle:
                table = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(table), 1)
            self.assertEqual(table[0]["source_attributes"], observation["source_attributes"])
            self.assertIn("</script><img", json.loads(observation["source_attributes"])[0]["attributes"]["Note"])
            html = (out / "comparison.html").read_text()
            self.assertNotIn("</script><img", html)
            payload = json.loads(re.search(r'<script id="gene-data" type="application/json">(.*?)</script>', html, re.S).group(1))
            self.assertEqual(payload["source_context"], gene["source_context"])
            self.assertIn("tables/gene_context.tsv", html)
            self.assertEqual(json.loads((out / "comparison.json").read_text()), summary)
            text = (out / "comparison.txt").read_text()
            self.assertIn("Source context: 1 of 1 unresolved CDS models", text)
            self.assertIn("gene-only pseudogene", text)
            self.assertIn("20 bp", text)
            self.assertIn("66.67% of 30 CDS bases", text)
            run_comparison(inputs, out)
            self.assertFalse((out / "tables/gene_context.tsv").exists())

    def test_no_context_exports_header_and_does_not_claim_absence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "genome.fna"
            fasta.write_text(">chr\n" + "ATG"*40 + "\n")
            inputs = []
            for label, kind in (("A", "CDS"), ("B", "region")):
                gff = root / f"{label}.gff3"
                gff.write_text(f"##gff-version 3\nchr\ttest\t{kind}\t10\t39\t.\t+\t0\tID={label}\n")
                inputs.append(GenomeInput(label, fasta, gff))
            out = root / "out"
            gene = run_comparison(inputs, out, gene_level=True)["gene_level_comparison"]
            self.assertEqual(gene["source_context"]["models_with_observations"], 0)
            self.assertEqual(len((out / "tables/gene_context.tsv").read_text().splitlines()), 1)
            self.assertIn("No qualifying source context recorded", (out / "comparison.txt").read_text())
