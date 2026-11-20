"""CDS extraction and translation utilities."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, MutableMapping
from dataclasses import replace

from Bio.Data.IUPACData import protein_letters_3to1

from annostat.models import CdsSequence, Feature


CALCULATION_POLICY = "annotation-metrics-1.6"
CDS_ANNOTATION_ATTRIBUTES = (
    "gene", "product", "Dbxref", "pseudo", "pseudogene", "partial",
    "start_range", "end_range", "transl_except", "exception",
)


def inherit_cds_attributes(features: tuple[Feature, ...]) -> tuple[Feature, ...]:
    """Use unambiguous parent context consistently, without editing inputs.

    A CDS inherits from a single gene/pseudogene parent, or from a single
    mRNA/transcript parent and that transcript's single gene/pseudogene parent.
    Precedence is CDS, then transcript, then gene. A gene's ``Name`` serves as
    its symbol only when it has no ``gene`` attribute and the name is not its ID
    or locus tag (funannotate convention). ``note`` is inherited so COG
    categories recorded there can be read.
    """
    genes = {feature.id: feature for feature in features if feature.type in {"gene", "pseudogene"}}
    transcripts = {feature.id: feature for feature in features if feature.type in {"mRNA", "transcript"}}

    def single_parent(feature: Feature, candidates: dict[str, Feature]) -> Feature | None:
        parent_ids = feature.attributes.get("Parent", "").split(",")
        return candidates.get(parent_ids[0]) if len(parent_ids) == 1 else None

    result = []
    for feature in features:
        inherited = {}
        if feature.type == "CDS":
            transcript = single_parent(feature, transcripts)
            gene = single_parent(transcript, genes) if transcript else single_parent(feature, genes)
            for parent in (gene, transcript):
                if parent is None:
                    continue
                inherited.update({key: parent.attributes[key] for key in (*CDS_ANNOTATION_ATTRIBUTES, "note")
                                  if key in parent.attributes})
                name = parent.attributes.get("Name", "").strip()
                if (parent is gene and name and "gene" not in parent.attributes
                        and name not in {parent.id, parent.attributes.get("locus_tag")}):
                    inherited["gene"] = name
                if parent.type == "pseudogene":
                    inherited["pseudo"] = "true"
        result.append(replace(feature, attributes={**inherited, **feature.attributes}) if inherited else feature)
    return tuple(result)


def complete_cds_start(feature: Feature) -> bool:
    """Whether initial phase and source attributes permit initiator recoding."""
    attributes = feature.attributes
    boundary = "end_range" if feature.strand == "-" else "start_range"
    generic_partial = attributes.get("partial", "").strip().lower() in {"1", "true", "yes"} and not ({"start_range", "end_range"} & attributes.keys())
    return not feature.phase and boundary not in attributes and not generic_partial


_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn")
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_GENETIC_CODE_OVERRIDES = {
    1: {},
    3: {"ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "TGA": "W"},
    4: {"TGA": "W"},
    11: {},
    25: {"TGA": "G"},
}
SUPPORTED_GENETIC_CODES = frozenset(_GENETIC_CODE_OVERRIDES)
BACTERIAL_GENETIC_CODES = frozenset({4, 11, 25})
_SUPPORTED_TEXT = "1, 3, 4, 11, or 25"
_START_CODONS = {
    1: {"TTG", "CTG", "ATG"},
    3: {"ATA", "ATG", "GTG"},
    4: {"TTA", "TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
    11: {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"},
    25: {"TTG", "ATG", "GTG"},
}
_TRANSL_EXCEPT_PATTERN = re.compile(
    r"pos:(?:([0-9]+)\.\.([0-9]+)|complement\(([0-9]+)\.\.([0-9]+)\)),aa:([A-Za-z*]+)",
    re.IGNORECASE,
)
_EXCEPTION_AMINO_ACIDS = {
    **{name.lower(): amino_acid for name, amino_acid in protein_letters_3to1.items()},
    "sec": "U",
    "pyl": "O",
    "term": "*",
    "other": "X",
    "*": "*",
}


def _declared_code(feature: Feature) -> int | None:
    """Return a feature's supported ``transl_table``; reject malformed or unsupported values."""

    if "transl_table" not in feature.attributes:
        return None
    raw_code = feature.attributes["transl_table"].strip()
    try:
        code = int(raw_code)
    except ValueError as error:
        raise ValueError(
            f"feature {feature.id!r} has invalid transl_table {raw_code!r}"
        ) from error
    if code not in SUPPORTED_GENETIC_CODES:
        raise ValueError(
            f"feature {feature.id!r} declares unsupported NCBI genetic code {code}; "
            f"choose {_SUPPORTED_TEXT}"
        )
    return code


def _cds_groups(features: Iterable[Feature]) -> tuple[dict[str, list[Feature]], dict[str, int]]:
    """Group CDS rows by explicit ID and collect per-sequence region tables."""

    grouped: dict[str, list[Feature]] = defaultdict(list)
    region_codes: dict[str, int] = {}
    for row_number, feature in enumerate(features):
        if feature.type == "region":
            code = _declared_code(feature)
            if code is not None:
                if region_codes.setdefault(feature.seqid, code) != code:
                    raise ValueError(
                        f"region rows for {feature.seqid!r} declare inconsistent transl_table values"
                    )
            continue
        if feature.type != "CDS":
            continue
        # Only explicit IDs join rows. Coordinate-based fallback IDs represent
        # independent anonymous features and must never merge accidentally.
        grouped[feature.attributes.get("ID") or f"__anonymous_{row_number}"].append(feature)
    return grouped, region_codes


def _cds_declared_code(group_id: str, parts: list[Feature], region_codes: Mapping[str, int]) -> int | None:
    """A CDS's own table, else its sequence region's table, else None."""

    codes = {code for code in map(_declared_code, parts) if code is not None}
    if len(codes) > 1:
        raise ValueError(f"multipart CDS {group_id!r} has inconsistent transl_table values")
    return next(iter(codes)) if codes else region_codes.get(parts[0].seqid)


def _table_declarations(features: Iterable[Feature]) -> tuple[frozenset[int], bool]:
    """Return declared tables and whether any CDS lacks a declaration."""

    features = list(features)
    grouped, region_codes = _cds_groups(features)
    declared = {code for feature in features if feature.type in {"CDS", "region"}
                for code in [_declared_code(feature)] if code is not None}
    undeclared = any(_cds_declared_code(group_id, parts, region_codes) is None
                     for group_id, parts in grouped.items())
    return frozenset(declared), undeclared


def resolve_genetic_code(
    features: Iterable[Feature],
    requested: int | None = None,
) -> int:
    """Resolve the default table for CDSs without their own ``transl_table``.

    A CDS's own ``transl_table`` (or its sequence region's) always applies to
    that CDS. The returned default covers the remaining CDSs: the explicit
    override, else a table declared uniformly by every CDS, else the bacterial
    default of 11. An override is rejected when every CDS declares a different
    table, because it would apply to nothing.
    """

    if requested is not None and requested not in SUPPORTED_GENETIC_CODES:
        raise ValueError(
            f"unsupported NCBI genetic code {requested}; choose {_SUPPORTED_TEXT}"
        )
    declared, undeclared = _table_declarations(features)
    if requested is not None:
        if declared and not undeclared and requested not in declared:
            raise ValueError(
                f"requested genetic code {requested} conflicts with GFF3 transl_table "
                f"{', '.join(map(str, sorted(declared)))}"
            )
        return requested
    if len(declared) == 1 and not undeclared:
        return next(iter(declared))
    return 11


def genetic_code_source(features: Iterable[Feature], requested: int | None = None) -> str:
    """Describe where the default table for undeclared CDSs came from."""

    if requested is not None:
        return "command_line"
    declared, undeclared = _table_declarations(features)
    return "gff3" if len(declared) == 1 and not undeclared else "default"


def recognized_start_codons(genetic_code: int) -> frozenset[str]:
    """Return the initiator codons defined for a supported translation table."""

    try:
        return frozenset(_START_CODONS[genetic_code])
    except KeyError as error:
        raise ValueError(
            f"unsupported NCBI genetic code {genetic_code}; choose {_SUPPORTED_TEXT}"
        ) from error


def reverse_complement(sequence: str) -> str:
    """Return an uppercase reverse complement, including IUPAC ambiguity codes."""

    return sequence.translate(_COMPLEMENT)[::-1].upper()


def extract_feature_sequence(feature: Feature, genome: Mapping[str, str]) -> str:
    """Extract a linear feature using GFF3's 1-based inclusive coordinates."""

    return extract_feature_sequence_with_topology(feature, genome, frozenset())


def extract_feature_sequence_with_topology(
    feature: Feature,
    genome: Mapping[str, str],
    circular_seqids: frozenset[str],
) -> str:
    """Extract and orient a feature, including circular-origin wrapping.

    Linear features extending beyond their reference sequence are rejected.
    Features on the minus strand are reverse-complemented after all referenced
    nucleotide segments have been joined.
    """

    try:
        reference = genome[feature.seqid]
    except KeyError as error:
        raise ValueError(
            f"GFF3 sequence {feature.seqid!r} is missing from the FASTA file"
        ) from error
    if feature.end > len(reference) and feature.seqid not in circular_seqids:
        raise ValueError(
            f"feature {feature.id!r} ends at {feature.end}, beyond "
            f"{feature.seqid!r} length {len(reference)}"
        )
    if feature.end <= len(reference):
        # Python's zero-based, end-exclusive slice needs only the start adjusted.
        sequence = reference[feature.start - 1 : feature.end]
    else:
        # Circular features may require more than one segment around the origin.
        remaining = feature.length
        position = (feature.start - 1) % len(reference)
        parts: list[str] = []
        while remaining:
            chunk_length = min(remaining, len(reference) - position)
            parts.append(reference[position : position + chunk_length])
            remaining -= chunk_length
            position = 0
        sequence = "".join(parts)
    return reverse_complement(sequence) if feature.strand == "-" else sequence.upper()


def translate_dna(
    sequence: str,
    *,
    codon_counts: Counter[str] | None = None,
    genetic_code: int = 11,
    complete_start: bool = True,
) -> str:
    """Translate complete codons using a supported NCBI code.

    Unknown or ambiguous codons become ``X``. Valid start codons of the selected
    table translate to methionine in the first position, and a terminal stop is
    omitted from the returned protein sequence. If supplied, ``codon_counts`` is
    updated during the same traversal used for translation.
    """

    if genetic_code not in _GENETIC_CODE_OVERRIDES:
        raise ValueError(f"unsupported NCBI genetic code {genetic_code}; choose {_SUPPORTED_TEXT}")
    sequence = sequence.upper()
    codon_table = {**_CODON_TABLE, **_GENETIC_CODE_OVERRIDES[genetic_code]}
    codons = [
        sequence[offset : offset + 3] for offset in range(0, len(sequence) - 2, 3)
    ]
    if codon_counts is not None:
        codon_counts.update(codons)
    amino_acids = [codon_table.get(codon, "X") for codon in codons]
    if complete_start and codons and codons[0] in _START_CODONS[genetic_code]:
        # NCBI translation tables define initiators separately from elongation.
        amino_acids[0] = "M"
    if amino_acids and amino_acids[-1] == "*":
        amino_acids.pop()
    return "".join(amino_acids)


def _apply_translation_exceptions(
    protein: str,
    parts: list[Feature],
    phase: int,
    sequence_length: int,
) -> tuple[str, tuple[int, ...], tuple[str, ...]]:
    """Apply position-specific NCBI ``transl_except`` annotations.

    The coordinate map is built in biological order, so the same logic handles
    plus, minus, multipart, and circular CDS features.
    """

    exception_text = ",".join(
        part.attributes.get("transl_except", "") for part in parts
        if part.attributes.get("transl_except")
    )
    if not exception_text:
        return protein, (), ()
    coordinates: list[int] = []
    for part in parts:
        raw_positions = (
            range(part.start, part.end + 1)
            if part.strand != "-"
            else range(part.end, part.start - 1, -1)
        )
        coordinates.extend((position - 1) % sequence_length + 1 for position in raw_positions)
    coding_coordinates = coordinates[phase:]
    amino_acids = list(protein)
    applied: list[int] = []
    unapplied: list[str] = []
    # Split only at a new exception, not at commas inside join/complement locations.
    entries = dict.fromkeys(re.split(r",\s*(?=\(?pos:)", exception_text, flags=re.IGNORECASE))
    for entry in entries:
        entry = entry.strip()
        normalized = re.sub(r"\s+", "", entry)
        if normalized.startswith("(") and normalized.endswith(")"):
            normalized = normalized[1:-1]
        match = _TRANSL_EXCEPT_PATTERN.fullmatch(normalized)
        if match is None:
            unapplied.append(entry)
            continue
        start, end = int(match.group(1) or match.group(3)), int(match.group(2) or match.group(4))
        replacement = _EXCEPTION_AMINO_ACIDS.get(match.group(5).lower())
        if (replacement is None or start < 1 or end - start != 2
                or end > max(sequence_length, max(part.end for part in parts))):
            unapplied.append(entry)
            continue
        handled = False
        target = {
            (position - 1) % sequence_length + 1
            for position in range(min(start, end), max(start, end) + 1)
        }
        for codon_index in range(len(coding_coordinates) // 3):
            codon_positions = set(coding_coordinates[codon_index * 3 : codon_index * 3 + 3])
            if codon_positions != target:
                continue
            if codon_index < len(amino_acids):
                amino_acids[codon_index] = replacement
                applied.append(codon_index)
                handled = True
                break
            # ``translate_dna`` removes a terminal stop. Restore it when an
            # explicit terminal exception recodes that codon as an amino acid.
            if codon_index == len(amino_acids):
                if replacement != "*":
                    amino_acids.append(replacement)
                    applied.append(codon_index)
                handled = True
                break
        if not handled:
            unapplied.append(entry)
    if amino_acids and amino_acids[-1] == "*" and len(amino_acids) == len(coding_coordinates) // 3:
        amino_acids.pop()
    return ("".join(amino_acids),
            tuple(sorted(index for index in set(applied) if index < len(amino_acids))),
            tuple(unapplied))


def extract_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
    *,
    genetic_code: int = 11,
) -> list[CdsSequence]:
    """Return extracted and translated records for every CDS feature."""

    return list(
        iter_cds_sequences(features, genome, circular_seqids, genetic_code=genetic_code)
    )


def iter_cds_sequences(
    features: Iterable[Feature],
    genome: Mapping[str, str],
    circular_seqids: frozenset[str] = frozenset(),
    *,
    codon_counts: Counter[str] | None = None,
    start_counts: MutableMapping[str, int] | None = None,
    genetic_code: int = 11,
) -> Iterable[CdsSequence]:
    """Yield CDS records while optionally accumulating codon statistics.

    The iterator keeps sequence export memory-efficient: callers can consume and
    write each CDS immediately instead of retaining every translated sequence.
    GFF3 phase is applied after strand orientation and before translation.
    Each CDS uses its own ``transl_table`` (or its sequence region's); others use
    ``genetic_code``.
    """

    grouped, region_codes = _cds_groups(features)
    for group_id, parts in grouped.items():
        code = _cds_declared_code(group_id, parts, region_codes) or genetic_code
        signatures = {(part.seqid, part.strand) for part in parts}
        if len(signatures) != 1:
            raise ValueError(
                f"multipart CDS {group_id!r} has inconsistent sequence IDs or strands"
            )
        ordered = sorted(
            parts,
            key=lambda item: item.start,
            reverse=parts[0].strand == "-",
        )
        feature = ordered[0]
        nucleotide = "".join(
            extract_feature_sequence_with_topology(part, genome, circular_seqids)
            for part in ordered
        )
        # Only the biologically first segment's phase trims the joined CDS.
        # Later phases describe frame continuity at segment boundaries.
        phase = feature.phase or 0
        coding_nucleotide = nucleotide[phase:]
        if start_counts is not None and len(coding_nucleotide) >= 3:
            start = coding_nucleotide[:3].upper()
            start_counts[start] = start_counts.get(start, 0) + 1
        protein = translate_dna(
            coding_nucleotide,
            codon_counts=codon_counts,
            genetic_code=code,
            complete_start=complete_cds_start(feature),
        )
        protein, exception_indices, unapplied = _apply_translation_exceptions(
            protein,
            ordered,
            phase,
            len(genome[feature.seqid]),
        )
        yield CdsSequence(
            feature=feature,
            nucleotide=nucleotide,
            coding_nucleotide=coding_nucleotide,
            protein=protein,
            segments=tuple(ordered),
            translation_exception_indices=exception_indices,
            unapplied_translation_exceptions=unapplied,
            genetic_code=code,
        )
