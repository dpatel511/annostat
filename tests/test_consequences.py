"""Measured protein consequences, distinct from CDS matching and biological truth."""
import unittest
from unittest.mock import patch

from annostat.consequences import protein_consequence


class ProteinConsequenceTests(unittest.TestCase):
    def test_identical_proteins_need_no_alignment(self):
        with patch("annostat.consequences.Align.PairwiseAligner", side_effect=AssertionError("unneeded alignment")):
            result = protein_consequence("MEQQ", "MEQQ")
        self.assertEqual(result["protein_consequence"], "identical")
        self.assertEqual(result["identical_aa"], 4)
        self.assertEqual(result["alignment_identity_percent"], 100)

    def test_gyrb_like_n_terminal_difference_has_exact_residue_evidence(self):
        result = protein_consequence("MAMEQQ", "MEQQ", "N")
        self.assertEqual(result["protein_consequence"], "terminal_only")
        self.assertEqual(result["n_terminal_delta_aa"], -2)
        self.assertEqual(result["c_terminal_delta_aa"], 0)
        self.assertEqual(result["n_terminal_residues_a"], "MA")
        self.assertEqual(result["n_terminal_residues_b"], "")
        self.assertEqual(result["identical_aa"], 4)
        self.assertEqual(result["substitution_count"], 0)
        self.assertEqual(result["gap_aa"], 2)
        self.assertAlmostEqual(result["alignment_identity_percent"], 100 * 4 / 6)
        self.assertEqual(result["alignment_b"], "--MEQQ")

    def test_c_terminal_extension_and_reversal(self):
        result = protein_consequence("MEQQ", "MEQQAC", "C")
        reverse = protein_consequence("MEQQAC", "MEQQ", "C")
        self.assertEqual(result["c_terminal_delta_aa"], 2)
        self.assertEqual(reverse["c_terminal_delta_aa"], -2)
        self.assertEqual(result["c_terminal_residues_b"], "AC")
        self.assertEqual(result["alignment_identity_percent"], reverse["alignment_identity_percent"])

    def test_coordinate_supported_end_resolves_repeated_terminal_sequence(self):
        n = protein_consequence("AAAA", "AAA", "N")
        c = protein_consequence("AAAA", "AAA", "C")
        self.assertEqual((n["n_terminal_delta_aa"], n["c_terminal_delta_aa"]), (-1, 0))
        self.assertEqual((c["n_terminal_delta_aa"], c["c_terminal_delta_aa"]), (0, -1))

    def test_internal_substitution_is_measured_not_called_functional_change(self):
        result = protein_consequence("MPEPTIDE", "MPEPAIDE")
        self.assertEqual(result["protein_consequence"], "translation_difference")
        self.assertEqual(result["substitution_count"], 1)
        self.assertEqual(result["identical_aa"], 7)
        self.assertEqual(result["gap_aa"], 0)
        self.assertEqual(result["alignment_identity_percent"], 87.5)
        self.assertIsNone(result["n_terminal_delta_aa"])

    def test_terminal_label_does_not_contradict_the_changed_genomic_boundary(self):
        # Initiator-M recoding can create a misleading prefix after a start shift.
        result = protein_consequence("MAA", "MA", "N")
        self.assertNotEqual(result["protein_consequence"], "terminal_only")

    def test_internal_gap_is_retained(self):
        result = protein_consequence("MPEPTIDE", "MPEPAAATIDE")
        self.assertEqual(result["protein_consequence"], "translation_difference")
        self.assertEqual(result["gap_aa"], 3)
        self.assertEqual(result["substitution_count"], 0)

    def test_unknown_empty_and_internal_stop_translations_are_limited(self):
        for left, right in (("MXX", "MXX"), ("", ""), ("MA*Q", "MAAQ"), (None, "MA")):
            with self.subTest(left=left, right=right):
                result = protein_consequence(left, right)
                self.assertEqual(result["protein_consequence"], "limited")
                self.assertIsNone(result["identical_aa"])
                self.assertIsNone(result["alignment_identity_percent"])

    def test_selenocysteine_can_be_compared_without_substitution_matrix_failure(self):
        result = protein_consequence("MUAQ", "MUAT")
        self.assertEqual(result["substitution_count"], 1)

    def test_multiple_optimal_placements_are_disclosed_not_silently_chosen(self):
        result = protein_consequence("AAAA", "AAA")
        self.assertEqual(result["protein_consequence"], "limited")
        self.assertTrue(result["alignment_ambiguous"])
        self.assertIsNone(result["substitution_count"])

    def test_large_complex_pairs_are_explicitly_limited_but_terminal_checks_scale(self):
        result = protein_consequence("M" * 2100 + "A", "M" * 2100 + "C")
        self.assertEqual(result["protein_consequence"], "limited")
        self.assertIn("4,000,000", result["consequence_note"])
        self.assertEqual(protein_consequence("MA" + "Q" * 10000, "Q" * 10000, "N")["identical_aa"], 10000)

    def test_known_codon_mapping_does_not_need_quadratic_alignment(self):
        keys = [(3*i+1, 3*i+2, 3*i+3) for i in range(2102)]
        with patch("annostat.consequences.Align.PairwiseAligner", side_effect=AssertionError("unneeded alignment")):
            result = protein_consequence("MA" + "Q"*2100, "Q"*2100, "N", codons_a=keys, codons_b=keys[2:])
        self.assertEqual(result["protein_consequence"], "terminal_only")
        self.assertEqual(result["shared_codon_count"], 2100)

    def test_repeated_or_reordered_codon_keys_are_not_arbitrarily_aligned(self):
        for left, right in (([(1,2,3)]*3, [(1,2,3)]*3), ([(1,2,3),(4,5,6),(7,8,9)], [(7,8,9),(4,5,6),(1,2,3)])):
            result = protein_consequence("MQQ", "MQQ", codons_a=left, codons_b=right)
            self.assertEqual(result["protein_consequence"], "limited")
            self.assertIsNone(result["alignment_identity_percent"])


if __name__ == "__main__":
    unittest.main()
