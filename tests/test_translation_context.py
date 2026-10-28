"""Independent translation oracles; genomic and protein positions are 1-based."""
from itertools import product
import random
import unittest

from Bio.Data import CodonTable, IUPACData
from Bio.Seq import Seq

from annostat.comparison import _jensen_shannon
from annostat.gene_comparison import cds_model, compare_cds
from annostat.models import Feature
from annostat.sequences import iter_cds_sequences, translate_dna


class TranslationContextTests(unittest.TestCase):
    def record(self, dna, *, phase=0, strand="+", attributes=None):
        genomic = str(Seq(dna).reverse_complement()) if strand == "-" else dna
        feature = Feature("chr", "test", "CDS", 1, len(dna), None, strand, phase,
                          {"ID": "cds", **(attributes or {})})
        return next(iter(iter_cds_sequences([feature], {"chr": genomic})))

    def test_all_canonical_codons_match_independent_ncbi_tables(self):
        for code in (4, 11, 25):
            table = CodonTable.unambiguous_dna_by_id[code]
            for bases in product("ACGT", repeat=3):
                codon = "".join(bases)
                dna = codon + "GCTTAA"
                for complete in (False, True):
                    with self.subTest(code=code, codon=codon, complete=complete):
                        expected = str(Seq(dna).translate(table=code)).removesuffix("*")
                        if complete and codon in table.start_codons:
                            expected = "M" + expected[1:]
                        self.assertEqual(translate_dna(dna, genetic_code=code,
                                                       complete_start=complete), expected)

    def test_nonzero_initial_phase_does_not_invent_an_initiator(self):
        for strand in ("+", "-"):
            for phase in (1, 2):
                with self.subTest(strand=strand, phase=phase):
                    record = self.record("A" * phase + "GTGAAATAA", phase=phase, strand=strand)
                    self.assertEqual(record.coding_nucleotide, "GTGAAATAA")
                    self.assertEqual(record.protein, "VK")

    def test_explicit_standard_amino_acid_exceptions_are_applied(self):
        for name, amino_acid in IUPACData.protein_letters_3to1.items():
            with self.subTest(amino_acid=name):
                record = self.record("ATGGGGTAA", attributes={"transl_except": f"(pos:4..6,aa:{name})"})
                self.assertEqual(record.protein, "M" + amino_acid)
                self.assertEqual(record.translation_exception_indices, (1,))

    def test_explicit_terminal_termination_is_not_an_internal_stop(self):
        record = self.record("ATGAAAGGG", attributes={"transl_except": "(pos:7..9,aa:TERM)"})
        self.assertEqual(record.protein, "MK")
        self.assertEqual(record.translation_exception_indices, ())

    def test_exception_handling_does_not_strip_a_penultimate_internal_stop(self):
        record = self.record("ATGAAATAATAG", attributes={"transl_except": "(pos:4..6,aa:Lys)"})
        self.assertEqual(record.protein, "MK*")

    def test_identical_first_codons_at_different_loci_explain_l_vs_m(self):
        dna = "TTGAAATTGAAATAA"
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                genome = {"chr": str(Seq(dna).reverse_complement()) if strand == "-" else dna}
                long = Feature("chr", "test", "CDS", 1, 15, None, strand, 0, {"ID": "long"})
                short = Feature("chr", "test", "CDS", 1 if strand == "-" else 7,
                                9 if strand == "-" else 15, None, strand, 0, {"ID": "short"})
                models = [cds_model(next(iter(iter_cds_sequences([f], genome))), 15, 11, False)
                          for f in (long, short)]
                row = compare_cds(models[:1], models[1:], ("A", "B"))["relationships"][0]
                self.assertEqual(row["protein_consequence"], "initiator_recoding")
                self.assertEqual((row["start_codon_a"], row["start_codon_b"]), ("TTG", "TTG"))
                self.assertEqual((row["first_amino_acid_a"], row["first_amino_acid_b"]), ("M", "M"))
                positions = "9,8,7" if strand == "-" else "7,8,9"
                self.assertEqual(row["first_codon_positions_b"], positions)
                self.assertIn(positions, row["consequence_note"])
                self.assertIn("A residue 3=L", row["consequence_note"])
                self.assertIn("B residue 1=M", row["consequence_note"])
                self.assertEqual(row["start_context_b"], "initiator")

    def test_first_codon_context_retains_partial_and_source_exception(self):
        for attributes, expected, amino_acid in (
            ({"partial": "true"}, "incomplete_start", "V"),
            ({"transl_except": "(pos:1..3,aa:Sec)"}, "translation_exception", "U"),
        ):
            with self.subTest(attributes=attributes):
                record = self.record("GTGAAATAA", attributes=attributes)
                model = cds_model(record, 9, 11, False)
                row = compare_cds([model], [model], ("A", "B"))["relationships"][0]
                self.assertEqual(row["start_context_a"], expected)
                self.assertEqual(row["first_amino_acid_a"], amino_acid)
                self.assertIsNone(model.initiator_amino_acid)

    def test_first_codon_positions_preserve_circular_origin_order(self):
        feature = Feature("chr", "test", "CDS", 8, 16, None, "+", 0, {"ID": "wrap"})
        record = next(iter(iter_cds_sequences([feature], {"chr": "GAAATAAAT"}, frozenset({"chr"}))))
        model = cds_model(record, 9, 11, True)
        row = compare_cds([model], [model], ("A", "B"))["relationships"][0]
        self.assertEqual(row["start_codon_a"], "ATG")
        self.assertEqual(row["first_codon_positions_a"], "8,9,1")

    def test_proportional_profiles_do_not_crash_jensen_shannon_distance(self):
        rng = random.Random(42)
        for _ in range(200):
            left = {str(j): rng.random() * 100 for j in range(3)}
            right = {key: value * 10 for key, value in left.items()}
            self.assertAlmostEqual(_jensen_shannon(left, right), 0, places=7)
        self.assertEqual(_jensen_shannon({"A": 100}, {"B": 100}), 1)


if __name__ == "__main__":
    unittest.main()
