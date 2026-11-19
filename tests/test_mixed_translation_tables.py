"""Per-CDS NCBI translation tables (GFF3 coordinates are 1-based inclusive).

RefSeq fungal GFF3 declares transl_table only on mitochondrial CDSs; nuclear CDSs
must use the selected or default table instead of inheriting the mitochondrial one.
"""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq

from annostat.comparison import GenomeInput, run_comparison
from annostat.models import Feature
from annostat.plots import write_bar_chart
from annostat.sequences import iter_cds_sequences, resolve_genetic_code, translate_dna
from annostat.validation import validate_annotation
from annostat.workflow import run_analysis

NUCLEAR = "ATGAAATGA"          # standard code: M K stop
MITO = "ATGTGAAAATAA"          # table 4: M W K stop; table 1 would give an internal stop


def cds(seqid, start, end, identifier, extra=None, phase=0):
    return Feature(seqid, "test", "CDS", start, end, None, "+", phase, {"ID": identifier, **(extra or {})})


class MixedTranslationTableTests(unittest.TestCase):
    def setUp(self):
        self.genome = {"chrI": NUCLEAR, "chrM": MITO}
        self.features = [cds("chrI", 1, 9, "nuclear"), cds("chrM", 1, 12, "mito", {"transl_table": "4"})]

    def test_declared_table_applies_only_to_its_cds(self):
        records = {r.feature.id: r for r in iter_cds_sequences(self.features, self.genome, genetic_code=1)}
        self.assertEqual((records["nuclear"].protein, records["nuclear"].genetic_code), ("MK", 1))
        self.assertEqual((records["mito"].protein, records["mito"].genetic_code), ("MWK", 4))
        self.assertEqual(records["nuclear"].protein, str(Seq(NUCLEAR).translate(table=1)).rstrip("*"))
        self.assertEqual(records["mito"].protein, str(Seq(MITO).translate(table=4)).rstrip("*"))

    def test_region_table_applies_to_cds_on_that_sequence(self):
        region = Feature("chrM", "test", "region", 1, 12, None, "+", None, {"ID": "chrM", "transl_table": "4"})
        features = [region, cds("chrI", 1, 9, "nuclear"), cds("chrM", 1, 12, "mito")]
        records = {r.feature.id: r for r in iter_cds_sequences(features, self.genome, genetic_code=1)}
        self.assertEqual((records["mito"].protein, records["nuclear"].protein), ("MWK", "MK"))

    def test_default_table_resolution(self):
        self.assertEqual(resolve_genetic_code(self.features, 1), 1)  # applies to undeclared CDSs
        self.assertEqual(resolve_genetic_code(self.features), 11)     # undeclared CDSs keep the default
        self.assertEqual(resolve_genetic_code(self.features[1:]), 4)  # uniform declaration
        with self.assertRaisesRegex(ValueError, "conflicts with GFF3"):
            resolve_genetic_code(self.features[1:], 11)               # override would apply to nothing

    def test_inconsistent_rows_of_one_cds_are_rejected(self):
        parts = [cds("chrM", 1, 6, "split", {"transl_table": "4"}),
                 cds("chrM", 7, 12, "split", {"transl_table": "11"})]
        with self.assertRaisesRegex(ValueError, "inconsistent transl_table"):
            list(iter_cds_sequences(parts, self.genome))

    def test_table_3_yeast_mitochondrial_code(self):
        self.assertEqual(translate_dna("ATACTGTGATAA", genetic_code=3), "MTW")
        self.assertEqual(translate_dna("ATACTGTGATAA", genetic_code=3, complete_start=False), "MTW")
        self.assertEqual(translate_dna("GTGCTGTAA", genetic_code=3), "MT")


class MixedTableWorkflowTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.fasta = self.root / "genome.fa"
        self.fasta.write_text(f">chrI\n{NUCLEAR}GTGAAATAA\n>chrM\n{MITO}\n", encoding="utf-8")
        self.gff = self.root / "genes.gff3"
        self.gff.write_text(
            "##gff-version 3\n"
            "chrI\ttest\tCDS\t1\t9\t.\t+\t0\tID=nuclear\n"
            "chrI\ttest\tCDS\t10\t18\t.\t+\t0\tID=gtg_start\n"
            "chrM\ttest\tCDS\t1\t12\t.\t+\t0\tID=mito;transl_table=4\n", encoding="utf-8")

    def test_validation_accepts_different_tables_on_different_cds(self):
        result = validate_annotation(self.fasta, self.gff)
        self.assertTrue(result["valid"])
        self.assertNotIn("TRANSLATION_TABLE", result["rule_counts"])

    def test_analysis_records_per_cds_tables_and_counts_starts_once(self):
        with patch("annostat.workflow.write_bar_chart", wraps=write_bar_chart) as chart:
            summary = run_analysis(self.fasta, self.gff, self.root / "out", "csv", genetic_code=1)
        self.assertEqual(summary["genetic_code"], 1)
        self.assertEqual(summary["genetic_code_counts"], {"1": 2, "4": 1})
        proteins = (self.root / "out" / "sequences" / "cds_protein.fasta").read_text(encoding="utf-8")
        self.assertIn(">mito", proteins)
        self.assertIn("MWK\n", proteins)
        self.assertNotIn("internal_stop_codon", summary["quality_control"]["issue_counts"])
        starts = next(call.args[2] for call in chart.call_args_list if call.args[1] == "First translated codons")
        # GTG is not a table-1 initiator: it is unrecognized, not also an ATG/GTG/TTG bar.
        self.assertEqual(starts["GTG"], 0)
        self.assertEqual(starts["ATG"], 2)
        self.assertEqual(starts["Unrecognized"], 1)
        self.assertEqual(sum(starts.values()), 3)
        self.assertEqual(summary["recognized_start_count"], 2)

    def test_gene_comparison_uses_each_models_table(self):
        dataset = lambda label: GenomeInput(label, self.fasta, self.gff, genetic_code=1)
        result = run_comparison([dataset("A"), dataset("B")], self.root / "cmp", gene_level=True)
        rows = {row["id_a"]: row for row in result["gene_level_comparison"]["relationships"]}
        self.assertEqual((rows["mito"]["genetic_code_a"], rows["nuclear"]["genetic_code_a"]), (4, 1))
        self.assertEqual(rows["mito"]["protein_consequence"], "identical")
        self.assertEqual(result["datasets"][0]["genetic_code_counts"], {"1": 2, "4": 1})


if __name__ == "__main__":
    unittest.main()
