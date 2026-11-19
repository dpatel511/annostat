"""Offline evidence-first reporting of accepted and unresolved CDS relationships."""
from __future__ import annotations

import json
from collections import Counter
from html import escape
from importlib.resources import files
from statistics import median

from annostat.gene_comparison import PRODUCT_STATUS_LABELS
from annostat.locus_context import CONTEXT_LABELS


CONSEQUENCE_LABELS = {
    "identical": "Identical protein",
    "terminal_only": "Terminal-only difference",
    "initiator_recoding": "Initiation-context difference",
    "translation_difference": "Other translated-sequence difference",
    "limited": "Limited protein evidence",
}

REVIEW_GROUPS = (
    ("Terminal boundary differences", {"terminal_only"},
     "Inspect the changed start/stop boundary and extra residues; neither model is automatically preferred."),
    ("Translation-context differences", {"initiator_recoding", "translation_difference"},
     "Inspect initiation recoding, translation tables, and source exceptions; these are not DNA mutations."),
    ("Limited protein evidence", {"limited"},
     "Inspect translation and codon-mapping limits, including unknown/internal-stop residues; unavailable counts are not zero."),
    ("Unresolved correspondence", {"unresolved"},
     "Review unmatched calls and candidate/ambiguous overlaps; these do not establish gene absence or split/merge events."),
)


def _review_groups(gene: dict):
    for title, categories, guidance in REVIEW_GROUPS:
        yield title, guidance, [row for row in gene["relationships"] if row["protein_consequence"] in categories]


RELATIONSHIP_LABELS = {
    "exact": "Exact CDS structure", "start_shift": "Alternative start boundary",
    "stop_shift": "Alternative stop boundary", "frame_disagreement": "Frame disagreement",
    "overlap": "Other overlap", "unmatched": "Unmatched call",
}


def review_loci(rows: list[dict]) -> list[dict]:
    """Group unresolved links sharing a CDS; accepted matches remain unchanged."""
    nodes, by_model = {}, {}
    for index, row in enumerate(rows):
        if row["status"] == "matched":
            continue
        keys = [(side, row["seqid"], row["strand"], row[f"id_{side}"])
                for side in ("a", "b") if row[f"id_{side}"] is not None]
        nodes[index] = keys
        for key in keys:
            by_model.setdefault(key, []).append(index)
    groups, seen = [], set()
    for index in nodes:
        if index in seen:
            continue
        pending, indices, models = [index], set(), set()
        while pending:
            current = pending.pop()
            if current in seen:
                continue
            seen.add(current)
            indices.add(current)
            for key in nodes[current]:
                if key not in models:
                    models.add(key)
                    pending.extend(by_model[key])
        kinds = {rows[i]["status"] for i in indices}
        kind = "ambiguous" if "ambiguous" in kinds else "candidate" if "candidate" in kinds else "unmatched"
        if kind == "ambiguous":
            reason = ("Repeated exact CDS structures prevent a unique ID assignment."
                      if all(rows[i]["relationship"] == "exact" for i in indices) else
                      "Multiple overlapping CDS calls remain; no unique one-to-one correspondence was accepted.")
        elif kind == "candidate":
            reason = "An overlap is retained, but it does not meet the accepted exact/boundary correspondence rules."
        else:
            reason = "No counterpart was retained for this CDS by the current correspondence rules."
        groups.append(dict(rows=sorted(indices), kind=kind, reason=reason,
                           cds_counts={side: sum(key[0] == side for key in models) for side in ("a", "b")}))
    return groups


def _description_percent(value: float | None) -> str:
    if value is None:
        return "Not available"
    return "<0.1%" if 0 < value < 0.1 else f"{value:.1f}%"


def _description_panel(gene: dict) -> str:
    description = gene["description_summary"]
    total = description["identical_protein_pairs"]
    labels = list(gene["cds_counts"])
    matrix = {(cell["product_status_a"], cell["product_status_b"]): cell
              for cell in description["product_status_matrix"]}
    rows = ""
    for a, title in PRODUCT_STATUS_LABELS.items():
        cells = ""
        for b in PRODUCT_STATUS_LABELS:
            cell = matrix[a, b]
            percentage = cell["percent_of_identical_pairs"]
            heat = min(1.0, (percentage or 0) / 100)
            cells += (f'<td><button type="button" class="gene-heat{" gene-heat-off" if a != b else ""}{" gene-heat-dark" if heat > .5 else ""}" '
                      f'style="--heat:{heat:.3f}" data-product-a="{a}" data-product-b="{b}" '
                      f'aria-label="A {title}; B {PRODUCT_STATUS_LABELS[b]}: {cell["pair_count"]:,} pairs" '
                      f'aria-controls="gene-table" {"disabled" if not cell["pair_count"] else ""}>'
                      f'<strong>{cell["pair_count"]:,}</strong><small>{escape(_description_percent(percentage))}</small></button></td>')
        rows += f'<tr><th scope="row">{title}</th>{cells}</tr>'
    changes = ", ".join(
        f'{cell["pair_count"]:,} {PRODUCT_STATUS_LABELS[cell["product_status_a"]].lower()} → {PRODUCT_STATUS_LABELS[cell["product_status_b"]].lower()}'
        for cell in description["product_status_matrix"]
        if cell["product_status_a"] != cell["product_status_b"] and cell["pair_count"])
    lead = (f'<p class="gene-lead"><b>{description["product_text_difference_pairs"]:,}</b> of {total:,} identical-protein pairs differ in product text'
            f'{f"; category changes: {changes}" if changes else ""}.</p>' if total else
            '<p class="gene-lead">No accepted pairs with reliable identical predicted proteins; label comparison not available.</p>')
    header_cells = "".join(f'<th scope="col">{title}</th>' for title in PRODUCT_STATUS_LABELS.values())
    return f'''<section id="gene-descriptions" class="gene-section" aria-labelledby="gene-labels-title">
<h2 id="gene-labels-title"><span>3</span>Product labels on identical proteins</h2>{lead}
<figure class="gene-figure-block"><div class="table-wrap"><table class="gene-product-matrix">
<thead><tr><th scope="col">A \\ B</th>{header_cells}</tr></thead><tbody>{rows}</tbody></table></div>
<figcaption><b>Table 3.</b> Product-label category in A ({escape(labels[0])}, rows) and B ({escape(labels[1])}, columns); percentages of {total:,} identical-protein pairs. Select a cell to list its pairs.</figcaption></figure></section>'''


PARTITION_LABELS = {
    "identical": "Paired, identical protein",
    "changed": "Paired, protein differs",
    "limited": "Paired, limited evidence",
    "unresolved": "Overlap, not paired",
    "unmatched": "No overlapping call",
}
_PARTITION_FILTERS = {
    "identical": 'data-filter-consequence="identical"',
    "changed": 'data-filter-consequence="changed"',
    "limited": 'data-filter-consequence="limited"',
    "unresolved": 'data-filter-status="candidate,ambiguous"',
    "unmatched": 'data-filter-structure="unmatched"',
}


def model_partition(gene: dict) -> list[tuple[str, Counter, int]]:
    """Assign every CDS model of each source to exactly one correspondence outcome."""
    result = []
    for side, label in zip("ab", gene["cds_counts"]):
        state: dict[tuple, str] = {}
        for row in gene["relationships"]:
            if row[f"id_{side}"] is None:
                continue
            key = (row["seqid"], row["strand"], row[f"id_{side}"])
            if row["status"] == "matched":
                consequence = row["protein_consequence"]
                state[key] = consequence if consequence in {"identical", "limited"} else "changed"
            else:
                state.setdefault(key, "unresolved" if row["status"] in {"candidate", "ambiguous"} else "unmatched")
        result.append((label, Counter(state.values()), gene["cds_counts"][label]))
    return result


def _share(count: int, total: int) -> str:
    if not total:
        return "Not available"
    value = 100 * count / total
    return "<0.1%" if 0 < value < 0.1 else f"{value:.1f}%"


def _partition_panel(partition: list[tuple[str, Counter, int]], matched: int) -> str:
    bars = ""
    for side, (label, counts, total) in zip("AB", partition):
        segments = "".join(
            f'<span class="gene-seg gene-seg-{key}" style="width:{100 * counts[key] / total:.3f}%" '
            f'title="{PARTITION_LABELS[key]}: {counts[key]:,} ({_share(counts[key], total)})"></span>'
            for key in PARTITION_LABELS if total and counts[key])
        description = "; ".join(f"{PARTITION_LABELS[key]} {counts[key]:,} ({_share(counts[key], total)})" for key in PARTITION_LABELS)
        bars += (f'<div class="gene-partition-row"><div class="gene-partition-label"><b>{side}</b> {escape(label)}'
                 f'<small>{matched:,} / {total:,} CDS models ({_share(matched, total)}) paired</small></div>'
                 f'<div class="gene-bar" role="img" aria-label="{escape(label)}: {escape(description)}">{segments}</div></div>')
    legend = ""
    for key, title in PARTITION_LABELS.items():
        values = "".join(f'<td>{counts[key]:,}</td><td class="gene-pct">{_share(counts[key], total)}</td>' for _, counts, total in partition)
        enabled = any(counts[key] for _, counts, _ in partition)
        legend += (f'<tr><th scope="row"><span class="gene-swatch gene-seg-{key}"></span>{title}</th>{values}'
                   f'<td><button type="button" class="gene-link" {_PARTITION_FILTERS[key]} aria-controls="gene-table" '
                   f'{"" if enabled else "disabled"}>View</button></td></tr>')
    names = "".join(f'<th scope="col" colspan="2">{side}<span class="gene-col-label"> · {escape(label)}</span></th>'
                    for side, (label, _, _) in zip("AB", partition))
    return (f'<figure class="gene-figure-block"><div class="gene-partition">{bars}</div>'
            f'<figcaption><b>Figure 1.</b> Outcome of every CDS model in each annotation; each model is counted once.</figcaption></figure>'
            f'<div class="table-wrap"><table class="gene-partition-table">'
            f'<thead><tr><th scope="col">Outcome</th>{names}<th scope="col"><span class="gene-sr">Rows</span></th></tr></thead>'
            f'<tbody>{legend}</tbody></table></div>')


_DELTA_BINS = ((None, -500, "≤−500"), (-499, -100, "−499 to −100"), (-99, -25, "−99 to −25"), (-24, -1, "−24 to −1"),
               (0, 0, "0"), (1, 24, "1 to 24"), (25, 99, "25 to 99"), (100, 499, "100 to 499"), (500, None, "≥500"))


def _length_histogram(deltas: list[int]) -> str:
    """Histogram of protein length change (B minus A, aa) as inline SVG."""
    counts = [sum((low is None or delta >= low) and (high is None or delta <= high) for delta in deltas)
              for low, high, _ in _DELTA_BINS]
    peak = max(counts) or 1
    width, height, top, base = 720, 210, 22, 166
    step = (width - 20) / len(counts)
    bars = ""
    for index, ((low, high, label), count) in enumerate(zip(_DELTA_BINS, counts)):
        x = 10 + index * step + 10
        bar_width = step - 20
        bar = (base - top) * count / peak
        kind = "same" if low == 0 and high == 0 else "shorter" if high is not None and high < 0 else "longer"
        centre = x + bar_width / 2
        bars += (f'<rect class="gene-hist-{kind}" x="{x:.1f}" y="{base - bar:.1f}" width="{bar_width:.1f}" height="{bar:.1f}">'
                 f'<title>{label} aa: {count:,}</title></rect>'
                 f'<text x="{centre:.1f}" y="{base - bar - 6:.1f}" text-anchor="middle" class="gene-hist-count">{count:,}</text>'
                 f'<text x="{centre:.1f}" y="{base + 16}" text-anchor="middle" class="gene-hist-label">{label}</text>')
    return (f'<svg class="gene-hist" viewBox="0 0 {width} {height}" role="img" aria-label="Protein length change, B minus A, '
            f'among {len(deltas):,} interpretable protein differences">{bars}'
            f'<line x1="10" x2="{width - 10}" y1="{base}" y2="{base}" class="gene-hist-axis"/>'
            f'<text x="10" y="{height - 4}" class="gene-hist-label">B shorter</text>'
            f'<text x="{width - 10}" y="{height - 4}" text-anchor="end" class="gene-hist-label">B longer</text>'
            f'<text x="{width / 2:.0f}" y="{height - 4}" text-anchor="middle" class="gene-hist-label">Δ length, aa (B − A)</text></svg>')


def _description_text(gene: dict) -> list[str]:
    description = gene["description_summary"]
    total = description["identical_protein_pairs"]
    lines = [f"Product-label comparison: {total:,} accepted pairs with reliable identical predicted proteins.",
             f'{description["product_text_difference_pairs"]:,} of {total:,} identical-protein pairs have different product text (normalized whitespace; missing and blank treated alike).']
    for cell in description["product_status_matrix"]:
        percentage = cell["percent_of_identical_pairs"]
        percent = _description_percent(percentage) if percentage is not None else "percentage unavailable"
        lines.append(f'A {PRODUCT_STATUS_LABELS[cell["product_status_a"]]} / B {PRODUCT_STATUS_LABELS[cell["product_status_b"]]}: '
                     f'{cell["pair_count"]:,} / {total:,} ({percent})')
    lines.append(f"Gene symbols among these same {total:,} pairs:")
    lines += [f'{label}: gene symbol provided {counts["provided"]:,} / {total:,}; missing {counts["missing"]:,}'
              for label, counts in description["gene_symbol_counts"].items()]
    lines.append("Unresolved CDS models (distinct): " + "; ".join(
        f"{label}: {count:,}" for label, count in description["unresolved_cds_counts"].items()))
    lines += [description["methods"], "Product-status counts: tables/gene_product_status.tsv"]
    return lines


def _jump(attribute: str) -> str:
    return f' <button type="button" class="gene-link" {attribute} aria-controls="gene-table">View rows</button>'


def render_gene_report(summary: dict) -> str:
    """Present scientific findings before the complete, offline evidence browser."""
    gene = summary["gene_level_comparison"]
    loci = review_loci(gene["relationships"])
    labels = list(gene["cds_counts"])
    matched = gene["status_counts"].get("matched", 0)
    counts = gene["consequence_counts"]
    identical = counts.get("identical", 0)
    changed = sum(counts.get(key, 0) for key in ("terminal_only", "initiator_recoding", "translation_difference"))
    deltas = [row["protein_length_delta_aa"] for row in gene["relationships"]
              if row["protein_consequence"] in {"terminal_only", "initiator_recoding", "translation_difference"}]
    limited = counts.get("limited", 0)
    structures = Counter(row["relationship"] for row in gene["relationships"] if row["status"] == "matched")
    status = gene["status_counts"]
    unresolved = len(gene["relationships"]) - matched
    agreement = gene["exact_agreement"]
    exact = f"{100 * agreement:.2f}%" if agreement is not None else "Not available"
    total_models = sum(gene["cds_counts"].values())
    partition = model_partition(gene)
    description = gene["description_summary"]
    ambiguous_groups = sum(locus["kind"] == "ambiguous" for locus in loci)

    def percent(count, total):
        return f"{100 * count / total:.1f}%" if total else "Not available"

    paired = " · ".join(f"{side} {percent(matched, total)}" for side, total in zip("AB", gene["cds_counts"].values()))
    unpaired = " · ".join(f"{side} {total - matched:,}" for side, total in zip("AB", gene["cds_counts"].values()))
    stats = "".join(f'<div><dt>{title}</dt><dd>{value}</dd><p>{note}</p></div>' for title, value, note in (
        ("Accepted pairs", f"{matched:,}", f"of CDS models paired: {paired}"),
        ("Identical proteins", percent(identical, matched),
         f"{identical:,} of {matched:,} accepted pairs ({percent(identical, matched)})" if matched else "No accepted pairs"),
        ("Structural agreement", exact, f'{gene["exact_model_count"]:,} shared CDS structures'),
        ("Models without a pair", unpaired, "unresolved overlap or no overlapping call"),
    ))
    if matched:
        lead = (f"{matched:,} CDS models pair one-to-one between the two annotations; "
                f"{percent(identical, matched)} of these pairs encode identical proteins.")
    else:
        lead = "No one-to-one CDS pairs were accepted."

    findings = [f'<b>{gene["exact_model_count"]:,}</b> CDS structures are identical in both annotations ({exact} structural agreement).'
                + (_jump('data-filter-structure="exact"') if structures["exact"] else "")]
    if structures["start_shift"] or structures["stop_shift"]:
        shift = (f'; B encodes the longer protein in <b>{sum(delta > 0 for delta in deltas):,}</b> and the shorter in '
                 f'<b>{sum(delta < 0 for delta in deltas):,}</b> (median |Δ| {median(abs(delta) for delta in deltas):g} aa)' if deltas else "")
        findings.append(f'<b>{structures["start_shift"]:,}</b> pairs differ at the start and <b>{structures["stop_shift"]:,}</b> at the stop{shift}.'
                        + (_jump('data-filter-consequence="changed"') if changed else ""))
    unmatched_counts = gene["unmatched_cds_counts"]
    if any(unmatched_counts.values()):
        findings.append(" and ".join(f'<b>{count:,}</b> {side}' for side, count in zip("AB", unmatched_counts.values()))
                        + " calls overlap no CDS in the other annotation." + _jump('data-filter-structure="unmatched"'))
    unresolved_models = [model_counts["unresolved"] for _, model_counts, _ in partition]
    if any(unresolved_models):
        findings.append(" and ".join(f'<b>{count:,}</b> {side}' for side, count in zip("AB", unresolved_models))
                        + f" models overlap but cannot be paired one-to-one ({ambiguous_groups:,} ambiguous groups)."
                        + _jump('data-filter-status="candidate,ambiguous"'))
    if description["product_text_difference_pairs"]:
        findings.append(f'<b>{description["product_text_difference_pairs"]:,}</b> identical proteins carry different product text.'
                        ' <a class="gene-link" href="#gene-descriptions">Table 3</a>')
    if limited:
        findings.append(f'<b>{limited:,}</b> pairs have limited protein evidence (internal stops, unknown residues or unmapped codons).'
                        + _jump('data-filter-consequence="limited"'))
    finding_items = "".join(f"<li>{item}</li>" for item in findings)

    outcomes = ""
    for key, label in CONSEQUENCE_LABELS.items():
        count = counts.get(key, 0)
        outcomes += (f'<tr data-consequence="{key}" aria-label="{label}: {count:,} of {matched:,} accepted pairs ({percent(count, matched)})">'
                     f'<th scope="row"><span class="gene-swatch gene-color-{key}"></span>{label}</th>'
                     f'<td class="gene-result-count"><strong>{count:,}</strong></td><td class="gene-pct">{percent(count, matched)}</td>'
                     f'<td class="gene-share"><span class="gene-mini gene-color-{key}" style="--share:{count / matched if matched else 0:.4f}"></span></td>'
                     f'<td><button type="button" class="gene-link" data-consequence="{key}" aria-controls="gene-table" {"disabled" if not count else ""}>View</button></td></tr>')
    boundary_rows = "".join(
        f'<tr><th scope="row"><span class="gene-swatch gene-color-{key}"></span>{label}</th><td>{structures[key]:,}</td>'
        f'<td class="gene-pct">{percent(structures[key], matched)}</td>'
        f'<td class="gene-share"><span class="gene-mini gene-color-{key}" style="--share:{structures[key] / matched if matched else 0:.4f}"></span></td></tr>'
        for key, label in (("exact", "Exact CDS structure"), ("start_shift", "Alternative start boundary"), ("stop_shift", "Alternative stop boundary")))
    length_block = '<p class="gene-muted">No interpretable protein differences.</p>'
    if deltas:
        length_block = f'''<figure class="gene-figure-block">{_length_histogram(deltas)}
<figcaption><b>Figure 2.</b> Protein length change among {len(deltas):,} interpretable differences. B shorter: <b>{sum(delta < 0 for delta in deltas):,}</b> · B longer: <b>{sum(delta > 0 for delta in deltas):,}</b> · Same length: <b>{deltas.count(0):,}</b>; median absolute Δ: <strong>{median(abs(delta) for delta in deltas):g} aa</strong>; signed range: {min(deltas):+d} to {max(deltas):+d} aa.</figcaption></figure>'''
    context = gene.get("source_context")
    context_line = (f' {context["models_with_observations"]:,} of {context["model_count"]:,} unresolved models have source-context observations '
                    '(<a href="tables/gene_context.tsv">TSV</a>).' if context is not None else "")
    payload = json.dumps(gene, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    locus_payload = json.dumps(loci, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    assets = files("annostat").joinpath("assets")
    css = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.css", "gene_report.css"))
    scripts = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.js", "gene_report.js"))
    license_text = escape(assets.joinpath("datatables-LICENSE.txt").read_text(encoding="utf-8"))
    consequence_options = "".join(f'<option value="{key}">{label}</option>' for key, label in CONSEQUENCE_LABELS.items())
    product_options = "".join(f'<option value="{a}|{b}">A {title_a} / B {title_b}</option>'
                              for a, title_a in PRODUCT_STATUS_LABELS.items() for b, title_b in PRODUCT_STATUS_LABELS.items())
    structure_options = "".join(f'<option value="{key}">{RELATIONSHIP_LABELS[key]}</option>'
                                for key in ("exact", "start_shift", "stop_shift", "frame_disagreement", "overlap", "unmatched"))
    definitions = "".join(f"<dt>{CONSEQUENCE_LABELS[key]}</dt><dd>{text}</dd>" for key, text in (
        ("identical", "Translated residues at corresponding genomic codons are identical."),
        ("terminal_only", "Only a terminal boundary differs; all shared residues match."),
        ("initiator_recoding", "Initiator-methionine recoding explains the shared-codon mismatch."),
        ("translation_difference", "Other translation or coding-path differences (genetic code, exceptions, partial annotations)."),
        ("limited", "Internal stops, unknown residues or unmapped codons prevent residue-level comparison."),
    ))
    return f'''<style>{css}</style><section class="gene-report" id="gene-report">
<nav class="gene-nav" aria-label="Report sections"><a href="#gene-overview">Summary</a><a href="#gene-differences">Differences</a><a href="#gene-descriptions">Product labels</a><a href="#gene-evidence">Evidence</a><a href="#gene-methods">Methods</a><a href="#report-data">Data</a></nav>
<section class="gene-section" id="gene-overview" aria-labelledby="gene-summary-title"><h2 id="gene-summary-title"><span>1</span>Summary</h2>
<p class="gene-lead">{lead}</p>
<dl class="gene-stats">{stats}</dl>
<div class="gene-two-col"><div>{_partition_panel(partition, matched)}</div>
<div><h3>Key findings</h3><ul class="gene-findings">{finding_items}</ul></div></div></section>
<section class="gene-section" id="gene-differences" aria-labelledby="gene-diff-title"><h2 id="gene-diff-title"><span>2</span>Differences between accepted pairs</h2>
<div class="gene-two-col gene-even"><div><div class="table-wrap"><table class="gene-outcomes"><caption><b>Table 1.</b> Protein comparison of {matched:,} accepted pairs.</caption>
<thead><tr><th scope="col">Outcome</th><th scope="col">Pairs</th><th scope="col">%</th><th scope="col"><span class="gene-sr">Share</span></th><th scope="col"><span class="gene-sr">Rows</span></th></tr></thead><tbody>{outcomes}</tbody></table></div></div>
<div><div class="table-wrap"><table class="gene-boundaries"><caption><b>Table 2.</b> CDS structure of accepted pairs.</caption>
<thead><tr><th scope="col">Relationship</th><th scope="col">Pairs</th><th scope="col">%</th><th scope="col"><span class="gene-sr">Share</span></th></tr></thead><tbody>{boundary_rows}</tbody></table></div></div></div>
<div class="gene-figure-wide">{length_block}</div>
<p class="gene-unresolved"><b>{unresolved:,} unresolved relationship rows</b>: {status.get("unmatched", 0):,} unmatched calls, {status.get("candidate", 0):,} candidate and {status.get("ambiguous", 0):,} ambiguous rows ({ambiguous_groups:,} ambiguous overlap groups).{context_line}
<button type="button" class="gene-link" data-consequence="unresolved" aria-controls="gene-table" {'disabled' if not unresolved else ''}>View rows</button></p></section>
{_description_panel(gene)}
<section class="gene-section" id="gene-evidence" aria-labelledby="gene-evidence-title"><h2 id="gene-evidence-title"><span>4</span>Evidence</h2>
<p class="gene-muted">All {len(gene["relationships"]):,} relationships. Select <b>Inspect</b> for coordinates, translations and alignment.</p>
<div class="gene-controls"><label>View <select id="gene-view"><option value="review">Differences / needs review</option><option value="all">All relationships</option></select></label><label>Protein evidence <select id="gene-consequence"><option value="">All relationships</option>{consequence_options}<option value="changed">Any interpretable protein difference</option><option value="unresolved">Unresolved correspondence</option></select></label>
<label>Product labels <select id="gene-product-status"><option value="">All label combinations</option>{product_options}</select></label>
<label>Structure <select id="gene-structure"><option value="">All structures</option>{structure_options}</select></label></div>
<p class="gene-muted" id="gene-view-summary" role="status" aria-live="polite"></p>
<noscript><p>JavaScript is required for the interactive table; the <a href="tables/gene_comparisons.tsv">TSV</a> and <a href="comparison.json">JSON</a> contain every relationship.</p></noscript>
<div class="table-wrap"><table id="gene-table" aria-label="Complete CDS relationships"></table></div>
<section id="gene-detail" class="gene-detail" aria-live="polite"><h3>Gene evidence</h3><p>Select a relationship to inspect its two annotation models.</p></section></section>
<section class="gene-section" id="gene-methods" aria-labelledby="gene-methods-title"><h2 id="gene-methods-title"><span>5</span>Methods</h2>
<dl class="gene-methods"><dt>Coordinates</dt><dd>GFF3 positions are 1-based inclusive. Deltas are B minus A: protein length in amino acids; boundary shifts in bp, positive downstream.</dd>
<dt>Pairing</dt><dd>CDS models are matched by structure, strand and frame, not names. Exact structures are paired first, then unique start/stop alternatives; other overlaps stay candidate or ambiguous.</dd>
<dt>Structural agreement</dt><dd>2 × {gene["exact_model_count"]:,} shared models / {total_models:,} CDS models across both sources.</dd>
{definitions}
<dt>Product labels</dt><dd>Missing (absent or blank), explicit hypothetical (contains “hypothetical protein”) or other text, compared only for identical proteins. Other text does not establish a known function.</dd>
<dt>Limitations</dt><dd>A descriptive comparison of two annotations. It does not identify the correct model, gene absence, mutations or functional effects.</dd></dl>
<details class="gene-fineprint"><summary>Full method notes</summary><p>{escape(gene["consequence_methods"])}</p><p>{escape(description["methods"])}</p><p>{escape(gene["scope"])}</p>{f"<p>{escape(context['methods'])}</p>" if context is not None else ""}</details>
<p class="downloads"><a href="comparison.txt">Summary (TXT)</a><a href="tables/gene_comparisons.tsv">All relationships (TSV)</a><a href="tables/gene_product_status.tsv">Product labels (TSV)</a><a href="plots/gene_consequences.svg">Consequence figure (SVG)</a></p>
<details class="gene-fineprint"><summary>Bundled table-library license</summary><p>DataTables 3.1.4, locally bundled; no network access.</p><pre class="gene-license">{license_text}</pre></details>
</section></section><script id="gene-data" type="application/json">{payload}</script><script id="gene-loci" type="application/json">{locus_payload}</script><script>{scripts}</script>'''


def gene_text_report(summary: dict) -> str:
    """Keep the human-readable summary short; the TSV retains every row."""
    gene = summary["gene_level_comparison"]
    loci = review_loci(gene["relationships"])
    matched = gene["status_counts"].get("matched", 0)
    agreement = gene["exact_agreement"]
    exact = f"{100 * agreement:.2f}%" if agreement is not None else "n/a (no CDS models)"
    labels = [profile["label"] for profile in summary["datasets"]]
    lines = [f'Annostat {summary["annostat_version"]} | CDS annotation comparison',
             f"A: {labels[0]}", f"B: {labels[1]}",
             "Scope: two annotations of identical DNA; genomic coordinates 1-based inclusive.",
             f"Accepted one-to-one CDS pairs: {matched:,}"]
    lines += [f"{label}: {gene['consequence_counts'].get(key, 0):,} / {matched:,} accepted pairs" for key, label in CONSEQUENCE_LABELS.items()]
    lines += [f"Exact structural agreement: {exact}; {gene['exact_model_count']:,} shared models (multiplicity included).",
              "Unmatched CDS calls: " + "; ".join(f"{label}: {count:,}" for label, count in gene["unmatched_cds_counts"].items()),
              f"Candidate relationship rows: {gene['status_counts'].get('candidate', 0):,}; ambiguous relationship rows: {gene['status_counts'].get('ambiguous', 0):,}.",
              "Relationship rows are not distinct genes. Unresolved pairs have no assigned protein consequence.", "",
              *_description_text(gene), "",
              "Review shortlist (up to 3 examples per group; all rows remain in TSV/JSON):",
              "Within each group, measured absolute protein-length differences are largest first; unavailable lengths follow. This is not a functional-impact ranking."]
    ambiguous = [locus for locus in loci if locus["kind"] == "ambiguous"]
    lines += ["", f"Ambiguous overlap groups: {len(ambiguous):,} (groups are not inferred genes)."]
    for locus in ambiguous[:3]:
        members = {side: sorted({gene["relationships"][i][f"id_{side}"] for i in locus["rows"]
                                 if gene["relationships"][i][f"id_{side}"] is not None}) for side in ("a", "b")}
        lines += [f"- A: {', '.join(members['a'])}; B: {', '.join(members['b'])}; "
                  f"{len(locus['rows'])} relationship rows. {locus['reason']}",
                  "  No automatic split/merge or pairwise protein consequence is inferred for these links."]
    if len(ambiguous) > 3:
        lines.append("Additional groups remain available in the HTML locus views and full relationship exports.")
    context = gene.get("source_context")
    if context is not None:
        lines += ["", f'Source context: {context["models_with_observations"]:,} of {context["model_count"]:,} unresolved CDS models have recorded observations.',
                  "Observations are source context, not causes or accepted CDS counterparts."]
        for observation in context["observations"][:5]:
            overlap = observation["overlap_bp"]
            measure = (f'; {overlap} bp physical overlap ({observation["cds_coverage_percent"]:.2f}% of '
                       f'{observation["cds_covered_bp"]} CDS bases; {observation["feature_coverage_percent"]:.2f}% of '
                       f'{observation["feature_covered_bp"]} source-feature bases)' if overlap is not None else
                       "; " + observation["note"])
            location = (f' at {observation["feature_segments"]} ({observation["feature_strand"]})'
                        if overlap is not None else "")
            lines.append(f'- {observation["dataset"]}: {observation["cds_id"]} | {CONTEXT_LABELS[observation["observation"]]} '
                         f'({observation["context_dataset"]}: {observation["feature_id"] or "sequence edge"}){location}{measure}.')
        if context["models_with_observations"] < context["model_count"]:
            lines.append("No qualifying source context recorded for some unresolved models; this does not establish biological absence.")
        if len(context["observations"]) > 5:
            lines.append("First 5 observations shown; all observations remain in tables/gene_context.tsv and JSON.")
        lines.append("Source-context methods: " + context["methods"])
    for title, guidance, rows in _review_groups(gene):
        lines += ["", f"{title}: {len(rows):,} relationship row{'s' if len(rows) != 1 else ''}", guidance]
        rows.sort(key=lambda row: (row["protein_length_delta_aa"] is None,
                                  -abs(row["protein_length_delta_aa"] or 0), row["seqid"],
                                  row["id_a"] or "", row["id_b"] or "", row["status"]))
        if not rows:
            lines.append("None in these annotations.")
        for row in rows[:3]:
            category = row["protein_consequence"]
            delta = row["protein_length_delta_aa"]
            measurement = f"B-A {delta:+d} aa" if delta is not None else "Protein length difference unavailable"
            if category not in {"limited", "unresolved"}:
                measurement += f"; {row['identical_aa']} aa identical; {row['substitution_count']} non-gap mismatches; {row['gap_aa']} gap aa"
            lines.append(f"- {row['gene_a'] or row['gene_b'] or 'unnamed'} | {row['id_a'] or 'no counterpart'} / {row['id_b'] or 'no counterpart'} | {row['seqid']} ({row['strand']}) | {row['status']} / {row['relationship']} | {CONSEQUENCE_LABELS.get(category, 'Unresolved correspondence')}; {measurement}.")
            lines.append(f"  A: {row['segments_a'] or 'no counterpart'}; B: {row['segments_b'] or 'no counterpart'} (1-based inclusive).")
            lines.append(f"  {row['consequence_note']}")
            for side in ("a", "b"):
                if row[f"id_{side}"] is not None:
                    lines.append(f"  {side.upper()} first codon: {row[f'start_codon_{side}'] or 'unavailable'}; "
                                 f"genomic bases {row[f'first_codon_positions_{side}'] or 'unavailable'}; "
                                 f"first predicted aa {row[f'first_amino_acid_{side}'] or 'unavailable'}; "
                                 f"context {row[f'start_context_{side}'] or 'unavailable'}.")
                    if row[f"unapplied_transl_except_{side}"]:
                        lines.append(f"  {side.upper()} unapplied transl_except: {row[f'unapplied_transl_except_{side}']}")
            for end in ("n", "c"):
                if row[f"{end}_terminal_delta_aa"]:
                    residues = row[f"{end}_terminal_residues_a"] or row[f"{end}_terminal_residues_b"]
                    lines.append(f"  {end.upper()} terminus B-A {row[f'{end}_terminal_delta_aa']:+d} aa; extra sequence in {'A' if row[f'{end}_terminal_residues_a'] else 'B'}: {residues[:60]}{'... (full sequence in TSV/JSON)' if len(residues) > 60 else ''}.")
    lines += ["", "Interpretation: differences do not establish functional effects or annotation correctness. Unmatched calls do not establish biological absence. Experimental support is not assessed by this report.",
              gene["consequence_methods"],
              "Full evidence: comparison.json; tables/gene_comparisons.tsv; tables/gene_product_status.tsv; tables/gene_context.tsv; interactive comparison.html; plots/gene_consequences.svg.", ""]
    return "\n".join(lines)
