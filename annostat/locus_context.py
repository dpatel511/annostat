"""Source observations for unresolved CDS models; never accept or correct matches."""
from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Sequence

from annostat.gene_comparison import CdsModel, _intervals
from annostat.models import Feature
from annostat.qc import _circular_intervals, _is_pseudogene


CONTEXT_COLUMNS = (
    "dataset", "seqid", "cds_id", "cds_strand", "cds_segments", "observation",
    "context_dataset", "feature_id", "feature_type", "feature_strand",
    "feature_segments", "feature_source", "overlap_bp", "cds_covered_bp",
    "feature_covered_bp", "cds_coverage_percent", "feature_coverage_percent",
    "source_attributes", "note",
)
CONTEXT_LABELS = {
    "gene_only_pseudogene_overlap": "Overlaps a source gene-only pseudogene record",
    "rna_overlap": "Overlaps a source RNA record",
    "opposite_strand_cds_overlap": "Overlaps an opposite-strand source CDS",
    "sequence_boundary": "Touches a sequence boundary",
}
CONTEXT_METHODS = (
    "Context is recorded for distinct unresolved CDS models, not relationship rows. "
    "Overlap observations use only the other annotation's records on the same sequence; "
    "RNA/pseudogene overlaps can have either or unspecified strand. Gene-only means no "
    "CDS descendant via Parent links. Coordinates are 1-based inclusive. Overlap and "
    "coverage denominators use unions of physical bases, including circular wrapping; "
    "CDS phase is not trimmed. Sequence-boundary context uses the CDS source's declared "
    "circularity. Source attributes retain decoded original rows, not inherited values. "
    "These observations do not establish the cause of disagreement, truncation, "
    "functional status, biological absence or a preferred annotation."
)


def _union(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _source_records(features: list[Feature], lengths: dict[str, int]) -> list[dict]:
    """Group explicit multipart IDs and identify genes with any CDS descendant."""
    parents = defaultdict(set)
    for feature in features:
        if feature.attributes.get("ID"):
            parents[feature.seqid, feature.id].update(feature.attributes.get("Parent", "").split(","))
    coding_ancestors = set()
    pending = [(feature.seqid, parent) for feature in features if feature.type == "CDS"
               for parent in feature.attributes.get("Parent", "").split(",") if parent]
    while pending:
        key = pending.pop()
        if key in coding_ancestors:
            continue
        coding_ancestors.add(key)
        pending.extend((key[0], parent) for parent in parents[key] if parent)
    grouped = defaultdict(list)
    for index, feature in enumerate(features):
        if feature.seqid not in lengths:
            continue
        kind = ("opposite_strand_cds_overlap" if feature.type == "CDS" else
                "rna_overlap" if feature.type.endswith("RNA") and feature.type != "mRNA" else
                "gene_only_pseudogene_overlap" if _is_pseudogene(feature)
                and (feature.type == "gene" or "pseudogene" in feature.type.lower())
                and (feature.seqid, feature.id) not in coding_ancestors else None)
        if kind:
            # Anonymous rows are separate observations, even at identical coordinates.
            key = (kind, feature.seqid, feature.strand, feature.type, feature.attributes.get("ID") or index)
            grouped[key].append(feature)
    records = []
    for (kind, seqid, strand, feature_type, _), parts in grouped.items():
        intervals = _union(interval for part in parts
                           for interval in _circular_intervals(part, lengths[seqid]))
        records.append(dict(kind=kind, seqid=seqid, strand=strand, id=parts[0].id,
                            type=feature_type, intervals=intervals,
                            source=";".join(sorted({part.source for part in parts})),
                            attributes=json.dumps([
                                dict(start=part.start, end=part.end, phase=part.phase,
                                     attributes=part.attributes) for part in parts
                            ], sort_keys=True)))
    return records


def collect_source_context(gene: dict, models: Sequence[list[CdsModel]],
                           features: Sequence[list[Feature]], labels: tuple[str, str]) -> dict:
    """Join unresolved models to source observations, without changing gene."""
    queries = []
    for side in (0, 1):
        suffix = "ab"[side]
        ids = {(row["seqid"], row["strand"], row[f"id_{suffix}"])
               for row in gene["relationships"] if row["status"] != "matched"
               and row[f"id_{suffix}"] is not None}
        queries.append([model for model in models[side]
                        if (model.feature.seqid, model.feature.strand, model.feature.id) in ids])
    lengths = {model.feature.seqid: model.sequence_length for source in models for model in source}
    records = [_source_records(source, lengths) for source in features]
    observations = []

    def observation(side, model, record, overlap, note):
        own = _intervals(model)
        total = sum(end-start+1 for start, end in own)
        other_total = sum(end-start+1 for start, end in record["intervals"]) if record else None
        return dict(dataset=labels[side], seqid=model.feature.seqid, cds_id=model.feature.id,
                    cds_strand=model.feature.strand,
                    cds_segments=";".join(f"{start}-{end}" for start, end, _ in model.blocks),
                    observation=record["kind"] if record else "sequence_boundary",
                    context_dataset=labels[1-side] if record else labels[side],
                    feature_id=record["id"] if record else None,
                    feature_type=record["type"] if record else None,
                    feature_strand=record["strand"] if record else None,
                    feature_segments=";".join(f"{start}-{end}" for start, end in record["intervals"]) if record else None,
                    feature_source=record["source"] if record else None,
                    overlap_bp=overlap, cds_covered_bp=total, feature_covered_bp=other_total,
                    cds_coverage_percent=100*overlap/total if overlap is not None else None,
                    feature_coverage_percent=100*overlap/other_total if overlap is not None else None,
                    source_attributes=record["attributes"] if record else None, note=note)

    for side in (0, 1):
        other = records[1-side]
        # Physical interval unions ensure each shared base contributes once.
        events = sorted(
            [(model.feature.seqid, start, end, 0, index)
             for index, model in enumerate(queries[side]) for start, end in _intervals(model)] +
            [(record["seqid"], start, end, 1, index)
             for index, record in enumerate(other) for start, end in record["intervals"]]
        )
        active, overlaps = ({}, {}), defaultdict(int)
        for seqid, start, end, category, index in events:
            for current in active:
                for key in list(current):
                    if key[0] != seqid or current[key] < start:
                        del current[key]
            for (_, other_index), other_end in active[1-category].items():
                q, r = (index, other_index) if category == 0 else (other_index, index)
                record, model = other[r], queries[side][q]
                if record["kind"] == "opposite_strand_cds_overlap" and (
                    {model.feature.strand, record["strand"]} != {"+", "-"}
                ):
                    continue
                overlaps[q, r] += min(end, other_end) - start + 1
            active[category][seqid, index] = end
        for (q, r), overlap in sorted(overlaps.items()):
            record = other[r]
            observations.append(observation(side, queries[side][q], record, overlap,
                                            CONTEXT_LABELS[record["kind"]] + "; context only, not an accepted CDS counterpart."))
        for model in queries[side]:
            intervals = _intervals(model)
            if intervals[0][0] != 1 and intervals[-1][1] != model.sequence_length:
                continue
            note = ("Source declares this sequence circular; the CDS touches its coordinate origin."
                    if model.circular else
                    "Circularity is not declared by this source; touching a sequence edge does not establish truncation.")
            observations.append(observation(side, model, None, None,
                                            note + f" Sequence length: {model.sequence_length} bp."))
    observations.sort(key=lambda row: (row["dataset"], row["seqid"], row["cds_strand"], row["cds_id"],
                                       row["observation"], row["feature_id"] or "", row["feature_segments"] or ""))
    return dict(model_count=sum(map(len, queries)),
                models_with_observations=len({(row["dataset"], row["seqid"], row["cds_strand"], row["cds_id"])
                                              for row in observations}),
                methods=CONTEXT_METHODS, observations=observations)
