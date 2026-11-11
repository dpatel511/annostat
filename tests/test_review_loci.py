"""Report grouping preserves unresolved graph relationships, not inferred genes."""
import unittest
from annostat.gene_report import review_loci


def row(a, b, status="ambiguous", relationship="overlap", seqid="chr", strand="+"):
    return dict(id_a=a, id_b=b, status=status, relationship=relationship, seqid=seqid, strand=strand)


class ReviewLocusTests(unittest.TestCase):
    def test_crossing_calls_form_one_group_with_unique_cds_counts(self):
        rows = [row("A", "x"), row("B", "x"), row("B", "y")]
        groups = review_loci(rows)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["rows"], [0, 1, 2])
        self.assertEqual(groups[0]["cds_counts"], {"a": 2, "b": 2})
        self.assertEqual(groups[0]["kind"], "ambiguous")
        self.assertNotIn("fusion", groups[0]["reason"])

    def test_many_to_one_is_not_two_independent_genes_or_two_loci(self):
        groups = review_loci([row("A", "x"), row("B", "x")])
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["cds_counts"], {"a": 2, "b": 1})

    def test_accepted_overlapping_pairs_are_not_merged_into_ambiguity(self):
        groups = review_loci([row("A", "x", "matched", "exact"), row("B", "y", "matched", "exact")])
        self.assertEqual(groups, [])

    def test_reused_ids_on_other_strands_or_contigs_do_not_join(self):
        groups = review_loci([row("A", "x"), row("A", "x", seqid="plasmid"), row("A", "x", strand="-")])
        self.assertEqual(len(groups), 3)

    def test_candidate_and_unmatched_calls_have_different_explanations(self):
        groups = review_loci([row("A", "x", "candidate", "frame_disagreement"), row("B", None, "unmatched", "unmatched")])
        self.assertEqual([g["kind"] for g in groups], ["candidate", "unmatched"])
        self.assertEqual(groups[1]["cds_counts"], {"a": 1, "b": 0})
        self.assertNotEqual(groups[0]["reason"], groups[1]["reason"])

    def test_duplicate_exact_structures_explain_id_ambiguity(self):
        group = review_loci([row("a1", "x", relationship="exact"), row("a2", "x", relationship="exact")])[0]
        self.assertIn("exact CDS structures", group["reason"])

    def test_empty_input_and_direction_swap(self):
        self.assertEqual(review_loci([]), [])
        group = review_loci([row("x", "A"), row("x", "B")])[0]
        self.assertEqual(group["cds_counts"], {"a": 1, "b": 2})
