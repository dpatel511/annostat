"""Spliced gene/mRNA/CDS models under NCBI table 1; GFF3 coordinates are 1-based inclusive.

Expected proteins come from Biopython ``CompoundLocation`` extraction with
0-based half-open coordinates, independently of Annostat's CDS assembly.
"""

from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from annostat.analysis import analyze_features
from annostat.comparison import GenomeInput, run_comparison
from annostat.parsers import parse_gff
from annostat.qc import feature_quality_findings
from annostat.sequences import inherit_cds_attributes, iter_cds_sequences

_STOPS = {"TAA", "TAG", "TGA"}
_SENSE = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in _STOPS]


def _no_internal_stop(coding: str) -> bool:
    return "*" not in str(Seq(coding[: len(coding) - len(coding) % 3]).translate(table=1))[:-1]


class FungalPanel:
    """Build a small genome with planted spliced models and their GFF3 rows."""

    def __init__(self) -> None:
        self.rng = random.Random(20261119)
        self.genome = list(self._filler(3000))
        self.rows: list[str] = []
        self.locations: dict[str, tuple[str, list[tuple[int, int]], int]] = {}

    def _filler(self, length: int) -> str:
        return "".join(self.rng.choice("ACGT") for _ in range(length))

    def _codons(self, count: int) -> str:
        return "".join(self.rng.choice(_SENSE) for _ in range(count))

    def place(self, start: int, sequence: str) -> None:
        self.genome[start - 1 : start - 1 + len(sequence)] = sequence

    def gene(self, gene_id: str, strand: str, transcripts: dict[str, list[tuple[int, int]]],
             *, first_phase: int = 0, product: str = "hypothetical protein", new_gene: bool = True) -> None:
        spans = [coordinate for blocks in transcripts.values() for block in blocks for coordinate in block]
        if new_gene:
            self.rows.append(f"chr\tfunannotate\tgene\t{min(spans)}\t{max(spans)}\t.\t{strand}\t.\tID={gene_id};")
        for transcript_id, blocks in transcripts.items():
            ordered = sorted(blocks, reverse=strand == "-")
            low = min(start for start, _ in blocks)
            high = max(end for _, end in blocks)
            self.rows.append(f"chr\tfunannotate\tmRNA\t{low}\t{high}\t.\t{strand}\t.\t"
                             f"ID={transcript_id};Parent={gene_id};product={product};")
            for number, (start, end) in enumerate(ordered, 1):
                self.rows.append(f"chr\tfunannotate\texon\t{start}\t{end}\t.\t{strand}\t.\t"
                                 f"ID={transcript_id}.exon{number};Parent={transcript_id};")
            coding_bases = None  # coding bases before this row, in transcription order
            for start, end in ordered:
                phase = first_phase if coding_bases is None else (3 - coding_bases % 3) % 3
                self.rows.append(f"chr\tfunannotate\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t"
                                 f"ID={transcript_id}.cds;Parent={transcript_id};")
                coding_bases = (end - start + 1 - first_phase if coding_bases is None
                                else coding_bases + end - start + 1)
            self.locations[f"{transcript_id}.cds"] = (strand, ordered, first_phase)

    def write(self, root: Path, label: str, rename: dict[str, str] | None = None) -> GenomeInput:
        fasta = root / "genome.fa"
        fasta.write_text(">chr\n" + "".join(self.genome) + "\n", encoding="utf-8")
        text = "##gff-version 3\n" + "\n".join(self.rows) + "\n"
        for old, new in (rename or {}).items():
            text = text.replace(old, new)
        gff = root / f"{label}.gff3"
        gff.write_text(text, encoding="utf-8")
        return GenomeInput(label, fasta, gff, genetic_code=1)

    def biopython_protein(self, cds_id: str) -> str:
        strand, ordered, first_phase = self.locations[cds_id]
        sign = 1 if strand == "+" else -1
        parts = [SimpleLocation(start - 1, end, strand=sign) for start, end in ordered]
        location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
        coding = str(SeqFeature(location).extract(Seq("".join(self.genome))))[first_phase:]
        # GenBank /translation may translate an unambiguous partial final codon;
        # Annostat excludes incomplete trailing bases, so compare whole codons.
        protein = str(Seq(coding[: len(coding) - len(coding) % 3]).translate(table=1))
        if first_phase == 0 and coding[:3] in {"ATG", "CTG", "TTG"}:
            protein = "M" + protein[1:]
        return protein.removesuffix("*")


def build_panel() -> FungalPanel:
    panel = FungalPanel()
    # Plus-strand gene: T1 uses exons 1-2-3, T2 skips exon 2 (60 bp, frame kept).
    exons = [(101, 161), (221, 280), (341, 402)]
    while True:
        coding = "ATG" + panel._codons(59) + "TAA"
        skipped = coding[:61] + coding[121:]
        if _no_internal_stop(skipped) and _no_internal_stop(coding):
            break
    for (start, end), offset in zip(exons, (0, 61, 121)):
        panel.place(start, coding[offset : offset + end - start + 1])
    panel.gene("MUM_000001", "+", {"MUM_000001-T1": exons, "MUM_000001-T2": [exons[0], exons[2]]})
    # Minus-strand gene with two exons; the biological first exon is at higher coordinates.
    coding = "CTG" + panel._codons(28) + "TAG"
    intron = "GT" + panel._filler(56) + "AG"
    panel.place(601, str(Seq(coding[:40] + intron + coding[40:]).reverse_complement()))
    panel.gene("MUM_000002", "-", {"MUM_000002-T1": [(711, 750), (601, 650)]})
    # 5'-partial model: initial phase 1, no initiator, and a 3' partial codon (2 trailing bases).
    coding = "C" + panel._codons(20) + "GC"
    panel.place(1001, coding)
    panel.gene("MUM_000003", "+", {"MUM_000003-T1": [(1001, 1000 + len(coding))]}, first_phase=1)
    # Internal TGA stop in a two-exon model: retained as '*', never removed.
    coding = "ATG" + panel._codons(10) + "TGA" + panel._codons(10) + "TAA"
    panel.place(1501, coding[:30])
    panel.place(1601, coding[30:])
    panel.gene("MUM_000004", "+", {"MUM_000004-T1": [(1501, 1530), (1601, 1600 + len(coding) - 30)]},
               product="kinase")
    return panel


RENAMED = {"MUM_000001": "FUN_9", "MUM_000002": "FUN_8", "MUM_000003": "FUN_7", "MUM_000004": "FUN_6"}


class FungalModelTests(unittest.TestCase):
    def setUp(self) -> None:
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.panel = build_panel()
        self.left = self.panel.write(self.root, "A")
        self.right = self.panel.write(self.root, "B", RENAMED)
        self.features = inherit_cds_attributes(tuple(parse_gff(self.left.gff)))

    def test_counts_separate_genes_transcripts_models_and_cds_rows(self) -> None:
        summary = analyze_features(self.features)
        self.assertEqual(summary["feature_type_counts"]["gene"], 4)
        self.assertEqual(summary["feature_type_counts"]["mRNA"], 5)
        self.assertEqual(summary["cds_count"], 5)
        self.assertEqual(summary["cds_row_count"], 3 + 2 + 2 + 1 + 2)
        self.assertEqual(summary["hypothetical_cds_count"], 4)

    def test_spliced_translations_match_biopython_on_both_strands(self) -> None:
        records = {record.feature.id: record for record in
                   iter_cds_sequences(self.features, {"chr": "".join(self.panel.genome)}, genetic_code=1)}
        self.assertEqual(set(records), set(self.panel.locations))
        for cds_id, record in records.items():
            with self.subTest(cds=cds_id):
                self.assertEqual(record.protein, self.panel.biopython_protein(cds_id))
        self.assertTrue(records["MUM_000002-T1.cds"].protein.startswith("M"))  # CTG initiator
        self.assertFalse(records["MUM_000003-T1.cds"].protein.startswith("M"))
        self.assertEqual(len(records["MUM_000003-T1.cds"].protein), 20)
        self.assertEqual(records["MUM_000004-T1.cds"].protein.count("*"), 1)
        self.assertEqual(records["MUM_000001-T1.cds"].feature.attributes["product"], "hypothetical protein")

    def test_isoforms_and_exon_rows_are_not_overlap_or_duplicate_findings(self) -> None:
        findings = feature_quality_findings(self.features, {"chr": "".join(self.panel.genome)})
        kinds = {finding.issue_type for finding in findings}
        self.assertFalse(kinds & {"contained_cds", "cds_overlap", "adjacent_duplicate_annotation"})

    def test_renamed_ids_with_unchanged_structures_match_completely(self) -> None:
        gene = run_comparison([self.left, self.right], self.root / "out", gene_level=True)["gene_level_comparison"]
        self.assertEqual(gene["cds_counts"], {"A": 5, "B": 5})
        self.assertEqual(gene["status_counts"], {"matched": 5})
        self.assertEqual(gene["exact_agreement"], 1.0)
        pairs = {row["id_a"]: row for row in gene["relationships"]}
        self.assertEqual(pairs["MUM_000001-T2.cds"]["id_b"], "FUN_9-T2.cds")
        self.assertEqual(pairs["MUM_000004-T1.cds"]["protein_consequence"], "limited")
        self.assertEqual(pairs["MUM_000001-T1.cds"]["protein_consequence"], "identical")
        self.assertEqual({row["genetic_code_a"] for row in gene["relationships"]}, {1})

    def test_cds_identical_isoforms_remain_ambiguous_not_forced_pairs(self) -> None:
        # A third isoform of the same gene with T1's CDS structure (e.g. differing only in UTRs).
        exons = [(101, 161), (221, 280), (341, 402)]
        self.panel.gene("MUM_000001", "+", {"MUM_000001-T3": exons}, new_gene=False)
        left = self.panel.write(self.root, "A")
        right = self.panel.write(self.root, "B", RENAMED)
        gene = run_comparison([left, right], self.root / "dup", gene_level=True)["gene_level_comparison"]
        ambiguous = {(row["id_a"], row["id_b"]) for row in gene["relationships"] if row["status"] == "ambiguous"}
        self.assertEqual(ambiguous, {(a, b) for a in ("MUM_000001-T1.cds", "MUM_000001-T3.cds")
                                     for b in ("FUN_9-T1.cds", "FUN_9-T3.cds")})
        self.assertEqual(gene["exact_model_count"], 6)
        self.assertEqual(gene["status_counts"]["matched"], 4)


class FungalReportTests(unittest.TestCase):
    """funannotate conventions: symbols on gene Name, COG categories in mRNA note."""

    def test_analysis_reads_funannotate_labels_and_skips_bacterial_rrna_expectations(self):
        from annostat.workflow import run_analysis
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel = build_panel()
            dataset = panel.write(root, "A")
            text = dataset.gff.read_text(encoding="utf-8")
            text = text.replace("\tgene\t101\t402\t.\t+\t.\tID=MUM_000001;", "\tgene\t101\t402\t.\t+\t.\tID=MUM_000001;Name=GLT1_1;")
            text = text.replace("ID=MUM_000002;", "ID=MUM_000002;Name=MUM_000002;")  # Name equal to ID is not a symbol
            text = text.replace("ID=MUM_000004-T1;Parent=MUM_000004;product=kinase;",
                                "ID=MUM_000004-T1;Parent=MUM_000004;product=kinase;note=EggNog:ENOG50,COG:S;")
            dataset.gff.write_text(text, encoding="utf-8")
            summary = run_analysis(dataset.fasta, dataset.gff, root / "out", "csv", genetic_code=1)
            report = (root / "out" / "report.html").read_text(encoding="utf-8")
        self.assertEqual(summary["cds_with_gene_count"], 2)  # both MUM_000001 isoforms
        self.assertEqual(summary["cds_with_cog_count"], 1)
        self.assertEqual(summary["cog_category_counts"], {"S": 1})
        self.assertEqual(summary["rna_counts"], {})  # mRNAs are not counted as RNA features
        self.assertNotIn("missing_rrna_type", summary["quality_control"]["issue_counts"])
        self.assertIn("Genome annotation report", report)
        self.assertNotIn("Bacterial annotation report", report)

    def test_yeast_style_trna_names_count_as_observed_amino_acids(self):
        from annostat.models import Feature
        from annostat.qc import _structural_rna_findings
        trnas = [Feature("chr", "t", "tRNA", 1, 72, None, "+", None, {"ID": f"rna{i}", "gene": name})
                 for i, name in enumerate(("tD(GUC)1", "tS(UGA)3"))]
        missing = {finding.details for finding in _structural_rna_findings(trnas, bacterial_rrna=False)}
        self.assertEqual(len(missing), 18)
        self.assertNotIn("no tRNA annotation found for Asp", missing)
        self.assertNotIn("no tRNA annotation found for Ser", missing)

    def test_gene_name_fallback_ignores_identifiers_and_cds_names(self):
        from annostat.models import Feature
        gene = Feature("chr", "t", "gene", 1, 9, None, "+", None, {"ID": "gene-b0001", "Name": "b0001", "locus_tag": "b0001"})
        cds = Feature("chr", "t", "CDS", 1, 9, None, "+", 0, {"ID": "cds1", "Parent": "gene-b0001", "Name": "hypothetical protein"})
        self.assertNotIn("gene", inherit_cds_attributes((gene, cds))[1].attributes)


if __name__ == "__main__":
    unittest.main()
