"""Offline evidence-first reporting of accepted and unresolved CDS relationships."""
from __future__ import annotations

import json
from html import escape
from importlib.resources import files
from statistics import median

from annostat.gene_comparison import PAIRED_OUTCOMES, PRODUCT_STATUS_LABELS
from annostat.locus_context import CONTEXT_LABELS


CHANGED_CONSEQUENCES = frozenset({"terminal_only", "initiator_recoding", "translation_difference"})
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
    if total:
        lead = (f'<p class="gene-lead"><b>{description["product_text_difference_pairs"]:,}</b> of {total:,} identical-protein pairs differ in product text'
                f'{f"; category changes: {changes}" if changes else ""}.</p>')
    else:
        lead = '<p class="gene-lead">No accepted pairs with reliable identical predicted proteins; label comparison not available.</p>'
    normalized = lambda value: " ".join((value or "").split())
    examples = [row for row in gene["relationships"] if row["status"] == "matched" and row["protein_consequence"] == "identical"
                and normalized(row["product_a"]) != normalized(row["product_b"])][:5]
    example_list = ""
    if examples:
        items = "".join(f'<tr><td>{escape(normalized(row["product_a"]) or "(missing)")}</td><td class="gene-arrow" aria-hidden="true">→</td>'
                        f'<td>{escape(normalized(row["product_b"]) or "(missing)")}</td><td><code>{escape(row["id_a"])}</code></td></tr>' for row in examples)
        example_list = (f'<div class="table-wrap"><table class="gene-examples"><caption>Examples of changed product text.</caption>'
                        f'<thead><tr><th scope="col">A</th><th scope="col"><span class="gene-sr">to</span></th><th scope="col">B</th><th scope="col">CDS (A)</th></tr></thead>'
                        f'<tbody>{items}</tbody></table></div>')
    header_cells = "".join(f'<th scope="col">{title}</th>' for title in PRODUCT_STATUS_LABELS.values())
    return f'''<section id="gene-descriptions" class="gene-section" aria-labelledby="gene-labels-title">
<h2 id="gene-labels-title"><span>4</span>Product labels of identical proteins</h2>{lead}{example_list}
<figure class="gene-figure-block"><div class="table-wrap"><table class="gene-product-matrix">
<thead><tr><th scope="col">A \\ B</th>{header_cells}</tr></thead><tbody>{rows}</tbody></table></div>
<figcaption><b>Table 5.</b> Product-label category in A ({escape(labels[0])}, rows) and B ({escape(labels[1])}, columns); percentages of {total:,} identical-protein pairs. Select a cell to list its pairs.</figcaption></figure></section>'''


def _share(count: int, total: int) -> str:
    if not total:
        return "Not available"
    value = 100 * count / total
    return "<0.1%" if 0 < value < 0.1 else f"{value:.1f}%"


def _outcome_table(gene: dict) -> str:
    outcomes = gene["model_outcomes"]
    labels = list(gene["cds_counts"])
    totals = list(gene["cds_counts"].values())
    counts = [outcomes["counts"][label] for label in labels]

    def row(key: str, title: str) -> str:
        values = [side[key] for side in counts]
        cells = "".join(f'<td class="gene-count-{side}"><strong>{value:,}</strong></td><td class="gene-pct">{_share(value, total)}</td>'
                        for side, value, total in zip("ab", values, totals))
        bars = "".join(f'<span class="gene-mini gene-mini-{side}" style="--share:{value / total if total else 0:.4f}"></span>'
                       for side, value, total in zip("ab", values, totals))
        return (f'<tr data-outcome="{key}"><th scope="row"><i class="gene-swatch gene-o-{key}"></i>{title}</th>{cells}<td class="gene-share">{bars}</td>'
                f'<td><button type="button" class="gene-link" data-filter-outcome="{key}" aria-controls="gene-table">View</button></td></tr>')

    body = ""
    for group, keys in (("Paired one-to-one", PAIRED_OUTCOMES),
                        ("Not paired", tuple(key for key in outcomes["labels"] if key not in PAIRED_OUTCOMES))):
        group_values = [sum(side[key] for key in keys) for side in counts]
        body += (f'<tr class="gene-group"><th scope="rowgroup">{group}</th>'
                 + "".join(f'<td>{value:,}</td><td class="gene-pct">{_share(value, total)}</td>' for value, total in zip(group_values, totals))
                 + '<td colspan="2"></td></tr>')
        body += "".join(row(key, outcomes["labels"][key]) for key in keys if any(side[key] for side in counts))
    head = "".join(f'<th scope="col" colspan="2">{side}<span class="gene-col-label"> · {escape(label)}</span></th>'
                   for side, label in zip("AB", labels))
    total_cells = "".join(f'<td>{total:,}</td><td></td>' for total in totals)
    return (f'<div class="table-wrap"><table class="gene-outcome-table"><caption><b>Table 1.</b> Outcome of every CDS model; '
            f'each model is counted once. Bars: A (upper), B (lower).</caption>'
            f'<thead><tr><th scope="col">Outcome</th>{head}<th scope="col"><span class="gene-sr">Share</span></th><th scope="col"><span class="gene-sr">Rows</span></th></tr></thead>'
            f'<tbody>{body}</tbody><tfoot><tr><th scope="row">CDS models</th>{total_cells}<td colspan="2"></td></tr></tfoot></table></div>')


def _profile_table(gene: dict) -> str:
    labels = list(gene["cds_counts"])
    profiles = gene["model_outcomes"]["profiles"]
    if not any(profiles[label]["no_overlap"]["models"] for label in labels):
        return '<p class="gene-muted">Every CDS model overlaps a CDS in the other annotation.</p>'
    groups = [(side, kind, profiles[label][kind]) for side, label in zip("AB", labels) for kind in ("no_overlap", "paired")]
    longest = max((item["median_protein_length_aa"] or 0) for _, _, item in groups) or 1

    def cell(value, text, share):
        if value is None:
            return '<td class="gene-muted">Not available</td>'
        return f'<td><span class="gene-value">{text}</span><span class="gene-meter" style="--share:{share:.4f}"></span></td>'

    head = "".join(f'<th scope="col">{side} {"no overlap" if kind == "no_overlap" else "paired"}</th>' for side, kind, _ in groups)
    rows = (f'<tr><th scope="row">CDS models</th>' + "".join(f'<td>{item["models"]:,}</td>' for _, _, item in groups) + "</tr>"
            + '<tr><th scope="row">Median protein length</th>' + "".join(
                cell(item["median_protein_length_aa"], f'{item["median_protein_length_aa"]:g} aa' if item["median_protein_length_aa"] is not None else "",
                     (item["median_protein_length_aa"] or 0) / longest) for _, _, item in groups) + "</tr>"
            + '<tr><th scope="row">Explicit “hypothetical protein”</th>' + "".join(
                cell(item["explicit_hypothetical_percent"], f'{item["explicit_hypothetical_percent"]:.1f}%' if item["explicit_hypothetical_percent"] is not None else "",
                     (item["explicit_hypothetical_percent"] or 0) / 100) for _, _, item in groups) + "</tr>")
    return (f'<div class="table-wrap"><table class="gene-profile"><caption><b>Table 4.</b> Models with no overlapping call in the other '
            f'annotation, compared with paired models of the same annotation.</caption><thead><tr><th scope="col"></th>{head}</tr></thead>'
            f'<tbody>{rows}</tbody></table></div>')


_FATE_SHORT = {
    "unchanged": "Unchanged", "start_changed": "Start changed", "stop_changed": "Stop changed",
    "duplicate_structure": "Identical, not unique", "boundary_alternative": "Boundary alt., not unique",
    "internal_change": "Internal structure", "partial_overlap": "Partial overlap",
    "frame_disagreement": "Frame disagreement", "no_overlap": "No overlap",
}


def _fate_figure(gene: dict) -> str:
    """Both annotations on one count scale; paired models are linked between the bars."""
    outcomes = gene["model_outcomes"]
    labels = list(gene["cds_counts"])
    totals = list(gene["cds_counts"].values())
    scale = max(totals) or 1
    keys = list(outcomes["labels"])

    def bar(side: str, label: str, total: int) -> str:
        counts = outcomes["counts"][label]
        segments = ""
        for key in keys:
            if not counts[key]:
                continue
            share = counts[key] / total
            text = f'<span>{_FATE_SHORT[key]} · {_share(counts[key], total)}</span>' if counts[key] / scale >= .14 else ""
            segments += (f'<i class="gene-o-{key}" style="width:{100 * share:.3f}%" '
                         f'title="{outcomes["labels"][key]}: {counts[key]:,} ({_share(counts[key], total)})">{text}</i>')
        return (f'<div class="gene-fate-row"><div class="gene-fate-name"><b class="gene-fate-side">{side}</b>'
                f'<span>{escape(label)}</span><small>{total:,} CDS models</small></div>'
                f'<div class="gene-fate-track" style="width:{100 * total / scale:.3f}%">{segments}</div></div>')

    paired = outcomes["counts"][labels[0]]
    links = "".join(f'<i class="gene-o-{key}" style="width:{100 * paired[key] / scale:.3f}%"></i>'
                    for key in PAIRED_OUTCOMES if paired[key])
    matched = sum(paired[key] for key in PAIRED_OUTCOMES)
    link_label = (f'<span class="gene-fate-link-label">{matched:,} models paired one-to-one</span>' if matched else "")
    legend = "".join(f'<li><i class="gene-o-{key}"></i>{outcomes["labels"][key]}</li>'
                     for key in keys if any(outcomes["counts"][label][key] for label in labels))
    description = "; ".join(f"{label}: " + ", ".join(f"{outcomes['labels'][key]} {outcomes['counts'][label][key]:,}"
                                                   for key in keys if outcomes["counts"][label][key]) for label in labels)
    return (f'<figure class="gene-fate" role="img" aria-label="{escape(description)}">'
            f'{bar("A", labels[0], totals[0])}'
            f'<div class="gene-fate-row gene-fate-links-row"><div class="gene-fate-name"></div><div class="gene-fate-links">{links}{link_label}</div></div>'
            f'{bar("B", labels[1], totals[1])}'
            f'<ul class="gene-fate-legend">{legend}</ul>'
            f'<figcaption><b>Figure 1.</b> Fate of every CDS model, drawn to one scale. Paired models (left) correspond one-to-one '
            f'between A and B; the remaining models overlap without a unique pair or have no overlapping call.</figcaption></figure>')


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


def render_gene_report(summary: dict) -> str:
    """Present each CDS model's outcome first, then supporting tables and the evidence browser."""
    gene = summary["gene_level_comparison"]
    loci = review_loci(gene["relationships"])
    labels = list(gene["cds_counts"])
    totals = list(gene["cds_counts"].values())
    outcome_counts = [gene["model_outcomes"]["counts"][label] for label in labels]
    matched = gene["status_counts"].get("matched", 0)
    counts = gene["consequence_counts"]
    deltas = [row["protein_length_delta_aa"] for row in gene["relationships"] if row["protein_consequence"] in CHANGED_CONSEQUENCES]
    status = gene["status_counts"]
    unresolved = len(gene["relationships"]) - matched
    agreement = gene["exact_agreement"]
    exact = f"{100 * agreement:.2f}%" if agreement is not None else "Not available"
    ambiguous_groups = sum(locus["kind"] == "ambiguous" for locus in loci)
    description = gene["description_summary"]

    def percent(count, total):
        return f"{100 * count / total:.1f}%" if total else "Not available"

    a, b = outcome_counts
    unpaired = [sum(side[key] for key in side if key not in PAIRED_OUTCOMES and key != "no_overlap") for side in outcome_counts]
    if totals[0]:
        lead = (f'Of the {totals[0]:,} CDS models in A, <b class="gene-t-unchanged">{_share(a["unchanged"], totals[0])}</b> are unchanged in B, '
                f'<b class="gene-t-start_changed">{_share(a["start_changed"] + a["stop_changed"], totals[0])}</b> differ only at the start or stop, '
                f'<b class="gene-t-partial_overlap">{_share(unpaired[0], totals[0])}</b> overlap a B model without a one-to-one pairing and '
                f'<b class="gene-t-no_overlap">{_share(a["no_overlap"], totals[0])}</b> have no overlapping call. '
                f'B has <b>{b["no_overlap"]:,}</b> calls with no overlapping A model.')
    else:
        lead = f'A has no CDS models; all {totals[1]:,} B models have no counterpart.'
    if not matched:
        lead = "No one-to-one CDS pairs were accepted. " + lead

    protein_rows = ""
    for key, label in CONSEQUENCE_LABELS.items():
        count = counts.get(key, 0)
        if not count:
            continue
        protein_rows += (f'<tr data-consequence="{key}" aria-label="{label}: {count:,} of {matched:,} accepted pairs ({percent(count, matched)})">'
                         f'<th scope="row">{label}</th><td class="gene-result-count"><strong>{count:,}</strong></td>'
                         f'<td class="gene-pct">{percent(count, matched)}</td>'
                         f'<td><button type="button" class="gene-link" data-consequence="{key}" aria-controls="gene-table">View</button></td></tr>')
    changed_pairs = a["start_changed"] + a["stop_changed"]
    unchanged_other = sum(row["status"] == "matched" and row["relationship"] == "exact" and row["protein_consequence"] != "identical"
                          for row in gene["relationships"])
    paired_lead = (f'<b>{changed_pairs:,}</b> of {matched:,} pairs ({percent(changed_pairs, matched)}) differ in structure: '
                   f'{a["start_changed"]:,} at the start and {a["stop_changed"]:,} at the stop.'
                   + (f' {unchanged_other:,} pairs with unchanged structure translate differently or have limited evidence.' if unchanged_other else
                      ' Every pair with unchanged structure encodes an identical protein.')) if matched else "No accepted pairs."
    direction_rows = ""
    for key, title in (("start_shift", "Start changed"), ("stop_shift", "Stop changed")):
        selected = [row["protein_length_delta_aa"] for row in gene["relationships"]
                    if row["status"] == "matched" and row["relationship"] == key and row["protein_consequence"] in CHANGED_CONSEQUENCES]
        if selected:
            direction_rows += (f'<tr><th scope="row">{title}</th><td>{sum(value > 0 for value in selected):,}</td>'
                               f'<td>{sum(value < 0 for value in selected):,}</td><td>{selected.count(0):,}</td></tr>')
    direction_table = (f'<div class="table-wrap"><table class="gene-direction"><caption><b>Table 3.</b> Protein length in B relative to A, among interpretable protein differences.</caption>'
                       f'<thead><tr><th scope="col">Pairs</th><th scope="col">B longer</th><th scope="col">B shorter</th><th scope="col">Same</th></tr></thead>'
                       f'<tbody>{direction_rows}</tbody></table></div>') if direction_rows else ""
    figure = ""
    if deltas:
        figure = f'''<figure class="gene-figure-block gene-figure-wide">{_length_histogram(deltas)}
<figcaption><b>Figure 2.</b> Protein length change among {len(deltas):,} interpretable differences. B shorter: <b>{sum(delta < 0 for delta in deltas):,}</b> · B longer: <b>{sum(delta > 0 for delta in deltas):,}</b> · Same length: <b>{deltas.count(0):,}</b>; median absolute Δ: <strong>{median(abs(delta) for delta in deltas):g} aa</strong>; signed range: {min(deltas):+d} to {max(deltas):+d} aa.</figcaption></figure>'''
    context = gene.get("source_context")
    context_line = (f' {context["models_with_observations"]:,} of {context["model_count"]:,} unpaired models have source-context '
                    'observations (overlapping RNA/pseudogene records, opposite-strand CDSs or sequence edges; <a href="tables/gene_context.tsv">TSV</a>).'
                    if context is not None and context["model_count"] else "")
    unpaired_lead = (f'<b>{totals[0] - matched:,}</b> A and <b>{totals[1] - matched:,}</b> B models have no one-to-one pair. '
                     'Overlapping but unpaired models are classified in Table 1 by their most similar relationship.' + context_line)
    payload = json.dumps(gene, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    locus_payload = json.dumps(loci, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    assets = files("annostat").joinpath("assets")
    css = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.css", "gene_report.css"))
    scripts = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.js", "gene_report.js"))
    license_text = escape(assets.joinpath("datatables-LICENSE.txt").read_text(encoding="utf-8"))
    consequence_options = "".join(f'<option value="{key}">{label}</option>' for key, label in CONSEQUENCE_LABELS.items())
    product_options = "".join(f'<option value="{x}|{y}">A {title_x} / B {title_y}</option>'
                              for x, title_x in PRODUCT_STATUS_LABELS.items() for y, title_y in PRODUCT_STATUS_LABELS.items())
    structure_options = "".join(f'<option value="{key}">{RELATIONSHIP_LABELS[key]}</option>'
                                for key in ("exact", "start_shift", "stop_shift", "frame_disagreement", "overlap", "unmatched"))
    outcome_definitions = "".join(f"<dt>{gene['model_outcomes']['labels'][key]}</dt><dd>{text}</dd>" for key, text in (
        ("unchanged", "Accepted pair with identical CDS structure (segments, strand and frame)."),
        ("start_changed", "Accepted pair sharing the stop and all other segments; the start differs."),
        ("stop_changed", "Accepted pair sharing the start and all other segments; the stop differs."),
        ("duplicate_structure", "Identical structure exists in the other annotation, but repeated structures (often isoforms) prevent a unique pairing."),
        ("boundary_alternative", "A start or stop alternative exists but is not unique."),
        ("internal_change", "Same start and stop as an overlapping model, different internal segments (e.g. intron positions)."),
        ("partial_overlap", "Frame-compatible overlap with different ends and/or segments."),
        ("frame_disagreement", "Overlap only in a different reading frame."),
        ("no_overlap", "No CDS on the same strand overlaps in the other annotation."),
    ))
    protein_definitions = "".join(f"<dt>{CONSEQUENCE_LABELS[key]}</dt><dd>{text}</dd>" for key, text in (
        ("identical", "Translated residues at corresponding genomic codons are identical."),
        ("terminal_only", "Only a terminal boundary differs; all shared residues match."),
        ("initiator_recoding", "Initiator-methionine recoding explains the shared-codon mismatch."),
        ("translation_difference", "Other translation or coding-path differences (genetic code, exceptions, partial annotations)."),
        ("limited", "Internal stops, unknown residues or unmapped codons prevent residue-level comparison."),
    ))
    return f'''<style>{css}</style><section class="gene-report" id="gene-report">
<nav class="gene-nav" aria-label="Report sections"><a href="#gene-overview">Summary</a><a href="#gene-differences">Paired models</a><a href="#gene-unpaired">Unpaired models</a><a href="#gene-descriptions">Product labels</a><a href="#gene-evidence">Evidence</a><a href="#gene-methods">Methods</a><a href="#report-data">Data</a></nav>
<section class="gene-section" id="gene-overview" aria-labelledby="gene-summary-title"><h2 id="gene-summary-title"><span>1</span>Summary</h2>
<p class="gene-lead">{lead}</p>{_fate_figure(gene)}{_outcome_table(gene)}</section>
<section class="gene-section" id="gene-differences" aria-labelledby="gene-diff-title"><h2 id="gene-diff-title"><span>2</span>Paired models that differ</h2>
<p class="gene-lead gene-lead-small">{paired_lead}</p>
<div class="gene-two-col gene-even"><div><div class="table-wrap"><table class="gene-outcomes"><caption><b>Table 2.</b> Protein comparison of {matched:,} accepted pairs.</caption>
<thead><tr><th scope="col">Protein</th><th scope="col">Pairs</th><th scope="col">%</th><th scope="col"><span class="gene-sr">Rows</span></th></tr></thead><tbody>{protein_rows}</tbody></table></div></div>
<div>{direction_table}</div></div>{figure}</section>
<section class="gene-section" id="gene-unpaired" aria-labelledby="gene-unpaired-title"><h2 id="gene-unpaired-title"><span>3</span>Unpaired models</h2>
<p class="gene-lead gene-lead-small">{unpaired_lead}</p>{_profile_table(gene)}
<p><button type="button" class="gene-link" data-consequence="unresolved" aria-controls="gene-table" {'disabled' if not unresolved else ''}>View all unpaired relationships</button></p></section>
{_description_panel(gene)}
<section class="gene-section" id="gene-evidence" aria-labelledby="gene-evidence-title"><h2 id="gene-evidence-title"><span>5</span>Evidence</h2>
<p class="gene-muted">All {len(gene["relationships"]):,} relationships. Select <b>Inspect</b> for coordinates, translations and alignment.</p>
<div class="gene-controls"><label>View <select id="gene-view"><option value="review">Differences / needs review</option><option value="all">All relationships</option></select></label><label>Protein evidence <select id="gene-consequence"><option value="">All relationships</option>{consequence_options}<option value="changed">Any interpretable protein difference</option><option value="unresolved">Unresolved correspondence</option></select></label>
<label>Product labels <select id="gene-product-status"><option value="">All label combinations</option>{product_options}</select></label>
<label>Structure <select id="gene-structure"><option value="">All structures</option>{structure_options}</select></label></div>
<p class="gene-muted" id="gene-view-summary" role="status" aria-live="polite"></p>
<noscript><p>JavaScript is required for the interactive table; the <a href="tables/gene_comparisons.tsv">TSV</a> and <a href="comparison.json">JSON</a> contain every relationship.</p></noscript>
<div class="table-wrap"><table id="gene-table" aria-label="Complete CDS relationships"></table></div>
<section id="gene-detail" class="gene-detail" aria-live="polite"><h3>Gene evidence</h3><p>Select a relationship to inspect its two annotation models.</p></section></section>
<section class="gene-section" id="gene-methods" aria-labelledby="gene-methods-title"><h2 id="gene-methods-title"><span>6</span>Methods</h2>
<dl class="gene-methods"><dt>Coordinates</dt><dd>GFF3 positions are 1-based inclusive. Deltas are B minus A: protein length in amino acids; boundary shifts in bp, positive downstream.</dd>
<dt>Pairing</dt><dd>CDS models are matched by structure, strand and frame, not names. Exact structures are paired first, then unique start/stop alternatives; other overlaps are not paired.</dd>
<dt>Model outcomes</dt><dd>Each model is counted once in its own annotation, by its most similar relationship (Table 1 order).</dd>{outcome_definitions}
{protein_definitions}
<dt>Product labels</dt><dd>Missing (absent or blank), explicit hypothetical (contains “hypothetical protein”) or other text, compared only for identical proteins. Other text does not establish a known function.</dd>
<dt>Relationship rows</dt><dd>{len(gene["relationships"]):,} rows: {matched:,} accepted, {status.get("unmatched", 0):,} unmatched, {status.get("candidate", 0):,} candidate and {status.get("ambiguous", 0):,} ambiguous ({ambiguous_groups:,} ambiguous overlap groups); {unresolved:,} unresolved relationship rows. One model can appear in several rows.</dd>
<dt>Structural agreement</dt><dd>{exact}: 2 × {gene["exact_model_count"]:,} shared models / {sum(totals):,} CDS models across both sources.</dd>
<dt>Limitations</dt><dd>A descriptive comparison of two annotations. It does not identify the correct model, gene gain or loss, mutations or functional effects.</dd></dl>
<details class="gene-fineprint"><summary>Full method notes</summary><p>{escape(gene["model_outcomes"]["methods"])}</p><p>{escape(gene["consequence_methods"])}</p><p>{escape(description["methods"])}</p><p>{escape(gene["scope"])}</p>{f"<p>{escape(context['methods'])}</p>" if context is not None else ""}</details>
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
    outcomes = gene.get("model_outcomes")
    if outcomes:
        lines.append("CDS model outcomes (each model counted once; A / B):")
        lines += [f"{title}: " + " / ".join(f"{counts[key]:,}" for counts in outcomes["counts"].values())
                  for key, title in outcomes["labels"].items()]
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
