"""CDS correspondence on identical assemblies, using 1-based inclusive intervals."""

from __future__ import annotations

import hashlib
from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import islice

from Bio.Data import CodonTable

from annostat.models import CdsSequence, Feature
from annostat.consequences import CONSEQUENCE_COLUMNS, protein_consequence
from annostat.sequences import CDS_ANNOTATION_ATTRIBUTES as _ATTRIBUTES, complete_cds_start, inherit_cds_attributes


GENE_COMPARISON_COLUMNS = (
    "dataset_a", "dataset_b", "seqid", "strand", "status", "relationship",
    *(f"{name}_{suffix}" for suffix in ("a", "b") for name in (
        "id", "segments", "phase", "cds_length", "protein_length", "genetic_code",
        *(attribute.lower() for attribute in _ATTRIBUTES),
    )),
    "overlap_bp", "start_delta_bp", "stop_delta_bp", "protein_identical", "changed_attributes",
    "start_codon_a", "start_codon_b", *CONSEQUENCE_COLUMNS,
)


def genome_fingerprint(genome: Mapping[str, str]) -> str:
    """Hash sequence IDs and uppercase bases independently of FASTA formatting."""

    digest = hashlib.sha256()
    for seqid in sorted(genome):
        sequence = genome[seqid]
        digest.update(seqid.encode("utf-8") + b"\0")
        digest.update(str(len(sequence)).encode("ascii") + b"\0")
        for offset in range(0, len(sequence), 65536):
            digest.update(sequence[offset:offset + 65536].upper().encode("ascii"))
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class CdsModel:
    """Retain coordinates, metadata, and translated protein evidence, not DNA.

    Each block stores its inclusive start/end and codon-position offset. Blocks
    remain in biological order; this also preserves repeated/slippage segments.
    """

    feature: Feature
    blocks: tuple[tuple[int, int, int], ...]
    length: int
    protein_length: int
    protein_sha256: str
    genetic_code: int
    circular: bool
    sequence_length: int
    protein: str | None = None
    start_codon: str | None = None
    initiator_amino_acid: str | None = None

    def protein_codons(self) -> list[tuple[int, int, int]]:
        """Map each predicted residue to its ordered 1-based genomic triplet."""
        positions = (position for start, end, _ in self.blocks for position in (range(start, end + 1) if self.feature.strand == "+" else range(end, start - 1, -1)))
        coding = islice(positions, self.feature.phase or 0, (self.feature.phase or 0) + 3 * self.protein_length)
        return list(zip(coding, coding, coding))

    @property
    def signature(self) -> tuple:
        """Identify a structural model independently of feature IDs and labels."""

        return (self.feature.seqid, self.feature.strand, self.feature.phase, self.blocks)

    @property
    def start(self) -> int:
        """Return the biological 5-prime boundary, not the lower coordinate."""

        return self.blocks[0][0 if self.feature.strand == "+" else 1]

    @property
    def stop(self) -> int:
        """Return the biological 3-prime boundary."""

        return self.blocks[-1][1 if self.feature.strand == "+" else 0]


def cds_model(record: CdsSequence, sequence_length: int, genetic_code: int, circular: bool) -> CdsModel:
    """Normalize circular coordinates and preserve the frame across CDS segments."""

    if record.feature.strand not in {"+", "-"}:
        raise ValueError("--gene-level requires explicit + or - strands for every CDS")
    direction = 1 if record.feature.strand == "+" else -1
    offset = -(record.feature.phase or 0)
    blocks = []
    for segment in record.segments or (record.feature,):
        intervals = [(segment.start, min(segment.end, sequence_length))]
        if circular and segment.end > sequence_length:
            intervals.append((1, segment.end - sequence_length))
        ordered_intervals = intervals if direction == 1 else reversed(intervals)
        for start, end in ordered_intervals:
            boundary = start if direction == 1 else end
            frame = (offset - direction * boundary) % 3
            adjacent = blocks and (blocks[-1][1] + 1 == start if direction == 1 else blocks[-1][0] - 1 == end)
            if adjacent and blocks[-1][2] == frame:
                previous = blocks.pop()
                blocks.append((min(previous[0], start), max(previous[1], end), frame))
            else:
                blocks.append((start, end, frame))
            offset += end - start + 1
    start_codon = record.coding_nucleotide[:3] if len(record.coding_nucleotide) >= 3 else None
    table = CodonTable.unambiguous_dna_by_id[genetic_code]
    initiator = table.forward_table.get(start_codon) if record.protein.startswith("M") and complete_cds_start(record.feature) and start_codon in table.start_codons else None
    return CdsModel(
        record.feature, tuple(blocks), record.length, len(record.protein),
        hashlib.sha256(record.protein.encode("ascii")).hexdigest(), genetic_code, circular,
        sequence_length, record.protein,
        start_codon, initiator,
    )


def _intervals(model: CdsModel) -> list[tuple[int, int]]:
    """Return disjoint physical intervals so overlap bases are not double-counted."""

    merged = []
    for start, end, _ in sorted(model.blocks):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _candidate_links(left: list[CdsModel], right: list[CdsModel]) -> tuple[list[set[int]], list[set[int]]]:
    """Sweep intervals to find same-sequence, same-strand cross-annotation overlaps."""

    links = ([set() for _ in left], [set() for _ in right])
    events = sorted(
        (model.feature.seqid, model.feature.strand, start, end, side, index)
        for side, models in enumerate((left, right))
        for index, model in enumerate(models)
        for start, end in _intervals(model)
    )
    active = [{}, {}]
    previous = None
    for seqid, strand, start, end, side, index in events:
        if previous != (seqid, strand):
            active = [{}, {}]
            previous = (seqid, strand)
        active = [{key: value for key, value in items.items() if value >= start} for items in active]
        for other in active[1 - side]:
            links[side][index].add(other)
            links[1 - side][other].add(index)
        active[side][index] = max(end, active[side].get(index, end))
    return links


def _boundary_alternative(short: CdsModel, long: CdsModel, at_start: bool) -> bool:
    """Require the shorter coding path to be a prefix or suffix of the longer."""

    if len(short.blocks) > len(long.blocks):
        return False
    shared = long.blocks[-len(short.blocks):] if at_start else long.blocks[:len(short.blocks)]
    if at_start:
        if short.blocks[1:] != shared[1:]:
            return False
        a, b = short.blocks[0], shared[0]
    else:
        if short.blocks[:-1] != shared[:-1]:
            return False
        a, b = short.blocks[-1], shared[-1]
    fixed_boundary = (1 if short.feature.strand == "+" else 0) if at_start else (0 if short.feature.strand == "+" else 1)
    return a[fixed_boundary] == b[fixed_boundary] and a[2] == b[2] and b[0] <= a[0] <= a[1] <= b[1]


def _relationship(left: CdsModel, right: CdsModel) -> str:
    """Classify structure without asserting that a candidate is the same gene."""

    if left.signature == right.signature:
        return "exact"
    shared_frames = [a[2] == b[2] for a in left.blocks for b in right.blocks if max(a[0], b[0]) <= min(a[1], b[1])]
    if not all(shared_frames):
        return "frame_disagreement"
    short, long = (left, right) if left.length <= right.length else (right, left)
    if left.stop == right.stop and _boundary_alternative(short, long, at_start=True):
        return "start_shift"
    if left.start == right.start and _boundary_alternative(short, long, at_start=False):
        return "stop_shift"
    return "overlap"


def _attribute(model: CdsModel, name: str) -> str | None:
    """Normalize whitespace and Dbxref ordering, retaining absent values as null."""

    value = model.feature.attributes.get(name)
    if value is None:
        return None
    if name == "Dbxref":
        return ",".join(sorted({item.strip() for item in value.split(",") if item.strip()}))
    return " ".join(value.split())


def _row(left: CdsModel | None, right: CdsModel | None, labels: tuple[str, str], status: str) -> dict[str, object]:
    """Build a flat, JSON/TSV-compatible relationship row."""

    model = left or right
    row = {"dataset_a": labels[0], "dataset_b": labels[1], "seqid": model.feature.seqid, "strand": model.feature.strand, "status": status}
    relationship = _relationship(left, right) if left and right else "unmatched"
    row["relationship"] = relationship
    for suffix, item in (("a", left), ("b", right)):
        row.update({
            f"id_{suffix}": item.feature.id if item else None,
            f"segments_{suffix}": ";".join(f"{start}-{end}" for start, end, _ in item.blocks) if item else None,
            f"phase_{suffix}": item.feature.phase if item else None,
            f"cds_length_{suffix}": item.length if item else None,
            f"protein_length_{suffix}": item.protein_length if item else None,
            f"genetic_code_{suffix}": item.genetic_code if item else None,
            **{f"{name.lower()}_{suffix}": _attribute(item, name) if item else None for name in _ATTRIBUTES},
        })
    overlap = None
    start_delta = stop_delta = None
    changed = None
    if left and right:
        overlap = sum(max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1) for a in _intervals(left) for b in _intervals(right))
        direction = 1 if left.feature.strand == "+" else -1
        if not (left.circular or right.circular):
            start_delta = direction * (right.start - left.start)
            stop_delta = direction * (right.stop - left.stop)
        elif relationship == "exact":
            start_delta = stop_delta = 0
        elif relationship in {"start_shift", "stop_shift"}:
            at_start = relationship == "start_shift"
            delta = direction * ((right.start - left.start) if at_start else (right.stop - left.stop))
            downstream = (right.length < left.length) if at_start else (right.length > left.length)
            delta = delta % left.sequence_length if downstream else -((-delta) % left.sequence_length)
            start_delta, stop_delta = (delta, 0) if at_start else (0, delta)
        differences = [name for name in _ATTRIBUTES if _attribute(left, name) != _attribute(right, name)]
        if left.genetic_code != right.genetic_code:
            differences.append("genetic_code")
        changed = ";".join(differences)
    row.update({
        "overlap_bp": overlap, "start_delta_bp": start_delta, "stop_delta_bp": stop_delta,
        "protein_identical": left.protein_sha256 == right.protein_sha256 if left and right else None,
        "changed_attributes": changed,
    })
    row.update(start_codon_a=left.start_codon if left else None, start_codon_b=right.start_codon if right else None)
    if status == "matched":
        boundary = {"start_shift": "N", "stop_shift": "C"}.get(relationship)
        row.update(protein_consequence(left.protein, right.protein, boundary,
                                       codons_a=left.protein_codons(), codons_b=right.protein_codons(),
                                       initiator_a=left.initiator_amino_acid, initiator_b=right.initiator_amino_acid))
    else:
        row.update(dict.fromkeys(CONSEQUENCE_COLUMNS))
        row.update(protein_consequence="unresolved", consequence_note="No accepted one-to-one CDS correspondence; no protein consequence assigned.")
    return row


def compare_cds(left: Iterable[CdsModel], right: Iterable[CdsModel], labels: tuple[str, str]) -> dict[str, object]:
    """Anchor exact/unique boundary matches; retain other candidates and ambiguity."""

    def sort_key(model: CdsModel) -> tuple:
        return (model.signature, model.feature.id, tuple(sorted(model.feature.attributes.items())))

    a, b = sorted(left, key=sort_key), sorted(right, key=sort_key)
    links_a, links_b = _candidate_links(a, b)
    signatures = (defaultdict(list), defaultdict(list))
    for side, models in enumerate((a, b)):
        for index, model in enumerate(models):
            signatures[side][model.signature].append(index)
    common = signatures[0].keys() & signatures[1].keys()
    exact_count = sum(min(len(signatures[0][key]), len(signatures[1][key])) for key in common)
    rows = []
    matched_a, matched_b = set(), set()

    def match(i: int, j: int) -> None:
        matched_a.add(i)
        matched_b.add(j)
        rows.append(_row(a[i], b[j], labels, "matched"))

    for key in sorted(common):
        if len(signatures[0][key]) == len(signatures[1][key]) == 1:
            match(signatures[0][key][0], signatures[1][key][0])
    # Treat start and stop alternatives together, so their priority cannot choose
    # between two different plausible counterparts.
    # Unresolved exact duplicates also block weaker boundary matches.
    eligible_a, eligible_b = defaultdict(set), defaultdict(set)
    for i, neighbours in enumerate(links_a):
        if i in matched_a or a[i].signature in common:
            continue
        for j in neighbours:
            if j not in matched_b and b[j].signature not in common and _relationship(a[i], b[j]) in {"start_shift", "stop_shift"}:
                eligible_a[i].add(j)
                eligible_b[j].add(i)
    for i, neighbours in sorted(eligible_a.items()):
        if len(neighbours) == 1:
            j = next(iter(neighbours))
            if len(eligible_b[j]) == 1:
                match(i, j)
    remaining_a = [neighbours - matched_b if i not in matched_a else set() for i, neighbours in enumerate(links_a)]
    remaining_b = [neighbours - matched_a if j not in matched_b else set() for j, neighbours in enumerate(links_b)]
    for i, neighbours in enumerate(remaining_a):
        if i in matched_a:
            continue
        if not neighbours:
            rows.append(_row(a[i], None, labels, "unmatched"))
        for j in sorted(neighbours):
            status = "candidate" if len(neighbours) == len(remaining_b[j]) == 1 else "ambiguous"
            rows.append(_row(a[i], b[j], labels, status))
    for j, neighbours in enumerate(remaining_b):
        if j not in matched_b and not neighbours:
            rows.append(_row(None, b[j], labels, "unmatched"))
    rows.sort(key=lambda row: (row["seqid"], row["strand"], row["id_a"] or "", row["id_b"] or "", row["status"]))
    return {
        "schema_version": "1.2", "coordinate_system": "1-based inclusive",
        "scope": "CDS models only; RNA and gene-only pseudogene records are not matched. Unmatched CDS calls do not establish biological absence.",
        "cds_counts": {labels[0]: len(a), labels[1]: len(b)},
        "exact_model_count": exact_count,
        "exact_agreement": 2 * exact_count / (len(a) + len(b)) if a or b else None,
        "status_counts": dict(sorted(Counter(row["status"] for row in rows).items())),
        "relationship_counts": dict(sorted(Counter(row["relationship"] for row in rows).items())),
        "consequence_counts": dict(sorted(Counter(row["protein_consequence"] for row in rows if row["status"] == "matched").items())),
        "unmatched_cds_counts": {label: sum(row["status"] == "unmatched" and row[f"id_{suffix}"] is not None for row in rows) for label, suffix in zip(labels, ("a", "b"))},
        "consequence_methods": "Accepted pairs use ordered genomic codon coordinates, not score-optimized protein alignment. Protein-length and terminal deltas are B minus A. Identity includes gap columns. Initiation-specific mismatches are separated from other non-gap mismatches; neither represents a DNA mutation on these identical assemblies. Unknown/stop residues and unresolved codon mapping are marked limited. Source partial/pseudogene/exception context is not experimental validation.",
        "relationships": rows,
    }
