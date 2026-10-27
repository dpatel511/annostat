"""Descriptive protein-sequence consequences; never functional-impact predictions."""
from __future__ import annotations

from Bio import Align


MAX_ALIGNMENT_CELLS = 4_000_000
ALIGNMENT_METHOD = "global; match 2; mismatch -1; gap open -5; gap extend -0.5"
CONSEQUENCE_COLUMNS = (
    "protein_consequence", "consequence_note", "protein_length_delta_aa",
    "identical_aa", "substitution_count", "gap_aa", "alignment_identity_percent",
    "n_terminal_delta_aa", "c_terminal_delta_aa",
    "n_terminal_residues_a", "n_terminal_residues_b",
    "c_terminal_residues_a", "c_terminal_residues_b",
    "alignment_method", "alignment_ambiguous", "alignment_a", "alignment_b",
    "shared_codon_count", "initiator_mismatch_count", "internal_mismatch_count",
)


def protein_consequence(a: str | None, b: str | None, boundary: str | None = None, *,
                        codons_a=None, codons_b=None, initiator_a=None, initiator_b=None) -> dict[str, object]:
    """Compare translations, preserving supplied genomic codon correspondence.

    Signed length changes are B minus A. Alignment identity includes gap columns.
    Without codon coordinates, only a descriptive string comparison is possible;
    multiple optimal fallback alignments are not assigned a residue-level verdict.
    """

    result = dict.fromkeys(CONSEQUENCE_COLUMNS)
    result.update(protein_consequence="limited", protein_length_delta_aa=len(b)-len(a) if a is not None and b is not None else None)
    if not a or not b:
        result["consequence_note"] = "Empty or unavailable translation."
        return result
    if set(a + b) - set("ACDEFGHIKLMNPQRSTVWYOU"):
        result["consequence_note"] = "Unknown residues or internal stop symbols; residue-level consequence not assigned."
        return result
    if codons_a is not None and codons_b is not None:
        return _coordinate_consequence(a, b, codons_a, codons_b, initiator_a, initiator_b, boundary, result)
    if a == b:
        result.update(protein_consequence="identical", identical_aa=len(a), substitution_count=0,
                      gap_aa=0, alignment_identity_percent=100.0, alignment_method="exact sequence equality",
                      n_terminal_delta_aa=0, c_terminal_delta_aa=0, consequence_note="Translated sequences are identical.")
        return result
    long, short = (a, b) if len(a) > len(b) else (b, a)
    ends = []
    if long.endswith(short):
        ends.append("N")
    if long.startswith(short):
        ends.append("C")
    end = boundary if boundary in ends else ends[0] if boundary is None and len(ends) == 1 else None
    if end:
        extra = len(long) - len(short)
        residues = long[:extra] if end == "N" else long[-extra:]
        aligned_short = "-" * extra + short if end == "N" else short + "-" * extra
        aligned_a, aligned_b = (long, aligned_short) if len(a) > len(b) else (aligned_short, long)
        result.update(protein_consequence="terminal_only", identical_aa=len(short), substitution_count=0,
                      gap_aa=extra, alignment_identity_percent=100 * len(short) / len(long),
                      n_terminal_delta_aa=len(b)-len(a) if end == "N" else 0,
                      c_terminal_delta_aa=len(b)-len(a) if end == "C" else 0,
                      n_terminal_residues_a=residues if end == "N" and len(a) > len(b) else "",
                      n_terminal_residues_b=residues if end == "N" and len(b) > len(a) else "",
                      c_terminal_residues_a=residues if end == "C" and len(a) > len(b) else "",
                      c_terminal_residues_b=residues if end == "C" and len(b) > len(a) else "",
                      alignment_method="exact suffix" if end == "N" else "exact prefix",
                      alignment_a=aligned_a, alignment_b=aligned_b,
                      consequence_note=f"{end}-terminal difference only; the shared {len(short):,} aa sequence is identical.")
        return result
    if len(a) * len(b) > MAX_ALIGNMENT_CELLS:
        result["consequence_note"] = "Complex alignment exceeds the 4,000,000-cell limit; no residue-level consequence assigned."
        return result
    aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                                   open_gap_score=-5, extend_gap_score=-0.5)
    alignments = iter(aligner.align(a, b))
    first = next(alignments)
    result["alignment_method"] = ALIGNMENT_METHOD
    if next(alignments, None) is not None:
        result.update(alignment_ambiguous=True, consequence_note="Multiple optimal alignments; gap placement and residue-level consequence remain unresolved.")
        return result
    left, right = str(first[0]), str(first[1])
    same = sum(x == y for x, y in zip(left, right))
    gaps = sum(x == "-" or y == "-" for x, y in zip(left, right))
    result.update(protein_consequence="translation_difference", identical_aa=same,
                  substitution_count=len(left)-same-gaps, gap_aa=gaps,
                  alignment_identity_percent=100 * same / len(left), alignment_ambiguous=False,
                  alignment_a=left, alignment_b=right,
                  consequence_note="Sequence differences are not explained by one exact terminal extension/shortening; inspect the alignment, not an inferred functional effect.")
    return result


def _coordinate_consequence(a, b, keys_a, keys_b, initiator_a, initiator_b, boundary, result):
    """Keep shared codons fixed; never optimize away an initiation mismatch."""
    result["alignment_method"] = "ordered genomic codon coordinates (1-based inclusive)"
    if len(keys_a) != len(a) or len(keys_b) != len(b) or len(set(keys_a)) != len(keys_a) or len(set(keys_b)) != len(keys_b):
        result["consequence_note"] = "Missing or repeated genomic codon correspondence; no residue-level consequence assigned."
        return result
    lookup = {key: index for index, key in enumerate(keys_b)}
    common = [(i, lookup[key]) for i, key in enumerate(keys_a) if key in lookup]
    if not common or any(right[1] <= left[1] for left, right in zip(common, common[1:])):
        result["consequence_note"] = "No unambiguous, ordered shared genomic codons; no residue-level consequence assigned."
        return result
    aligned_a, aligned_b = [], []
    previous_a = previous_b = 0
    initiation = mismatches = 0
    for i, j in common + [(len(a), len(b))]:
        extra_a, extra_b = a[previous_a:i], b[previous_b:j]
        if extra_a and extra_b:
            result["consequence_note"] = "Unshared codons on both sides of a mapping interval; no arbitrary residue correspondence assigned."
            return result
        aligned_a.extend((extra_a, "-" * len(extra_b)))
        aligned_b.extend(("-" * len(extra_a), extra_b))
        if i < len(a):
            aligned_a.append(a[i]); aligned_b.append(b[j])
            if a[i] != b[j]:
                is_initiation = (i == 0 and a[i] == "M" and b[j] == initiator_a) or (j == 0 and b[j] == "M" and a[i] == initiator_b)
                initiation += bool(is_initiation)
                mismatches += not is_initiation
        previous_a, previous_b = i + 1, j + 1
    left, right = "".join(aligned_a), "".join(aligned_b)
    gaps = len(left) - len(common)
    same = len(common) - initiation - mismatches
    first_a, first_b = common[0]; last_a, last_b = common[-1]
    n_delta = first_b - first_a
    c_delta = (len(b)-last_b-1) - (len(a)-last_a-1)
    internal_gaps = gaps - first_a - first_b - (len(a)-last_a-1) - (len(b)-last_b-1)
    boundary_explains_gaps = (boundary == "N" and c_delta == 0) or (boundary == "C" and n_delta == 0)
    category = "translation_difference" if mismatches or internal_gaps or (gaps and not boundary_explains_gaps) else "initiator_recoding" if initiation else "terminal_only" if gaps else "identical"
    notes = {
        "identical": "Translated sequences at corresponding genomic codons are identical.",
        "terminal_only": f"Terminal boundary difference only; all {same:,} shared translated residues are identical.",
        "initiator_recoding": f"Boundary/initiation-context difference; {initiation} mismatch(es) follow initiator-methionine recoding. The other {same:,} shared residues are identical, not evidence of internal DNA mutations.",
        "translation_difference": "Other differences at corresponding genomic codons or coding paths; inspect translation tables, exceptions, and partial annotations. These identical assemblies do not establish DNA mutations or functional effects.",
    }
    result.update(protein_consequence=category, consequence_note=notes[category],
                  identical_aa=same, substitution_count=initiation+mismatches, gap_aa=gaps,
                  alignment_identity_percent=100*same/len(left), alignment_ambiguous=False,
                  shared_codon_count=len(common), initiator_mismatch_count=initiation, internal_mismatch_count=mismatches,
                  n_terminal_delta_aa=n_delta, c_terminal_delta_aa=c_delta,
                  n_terminal_residues_a=a[:first_a], n_terminal_residues_b=b[:first_b],
                  c_terminal_residues_a=a[last_a+1:], c_terminal_residues_b=b[last_b+1:],
                  alignment_a=left if category != "identical" else None,
                  alignment_b=right if category != "identical" else None)
    return result
