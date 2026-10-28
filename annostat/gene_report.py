"""Offline evidence-first reporting of accepted and unresolved CDS relationships."""
from __future__ import annotations

import json
from collections import Counter
from html import escape
from importlib.resources import files
from pathlib import Path
from statistics import median


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


def render_gene_report(summary: dict, plot: Path) -> str:
    """Present scientific findings before the complete, offline evidence browser."""
    gene = summary["gene_level_comparison"]
    matched = gene["status_counts"].get("matched", 0)
    counts = gene["consequence_counts"]
    identical = counts.get("identical", 0)
    changed = sum(counts.get(key, 0) for key in ("terminal_only", "initiator_recoding", "translation_difference"))
    deltas = [row["protein_length_delta_aa"] for row in gene["relationships"]
              if row["protein_consequence"] in {"terminal_only", "initiator_recoding", "translation_difference"}]
    limited = counts.get("limited", 0)
    structures = Counter(row["relationship"] for row in gene["relationships"] if row["status"] == "matched")
    unresolved = len(gene["relationships"]) - matched
    agreement = gene["exact_agreement"]
    exact = f"{100 * agreement:.2f}%" if agreement is not None else "Not available"

    def percent(count, total):
        return f"{100 * count / total:.1f}%" if total else "Not available"

    cards = ""
    for title, count in (("Identical proteins", identical), ("Interpretable protein differences", changed), ("Limited protein evidence", limited)):
        denominator = f"{count:,} of {matched:,} accepted pairs ({percent(count, matched)})" if matched else "No accepted pairs; percentage unavailable"
        percentage = f"{percent(count, matched)} of accepted pairs" if matched else "Percentage not available"
        cards += f'<div aria-label="{denominator}"><span>{title}</span><strong>{count:,}</strong><small>{percentage}</small></div>'
    sources = ""
    for side, (label, total) in zip(("A", "B"), gene["cds_counts"].items()):
        sources += f'''<article class="gene-source"><span class="gene-source-label">{side} · {escape(label)}</span>
<strong>{total:,} <small>CDS models</small></strong>
<p>Accepted counterpart: <b>{matched:,} / {total:,} CDS models ({percent(matched, total)})</b></p>
<p class="note">{total - matched:,} models are not resolved one-to-one, including {gene["unmatched_cds_counts"][label]:,} unmatched calls.</p></article>'''
    if not matched:
        headline = "No one-to-one CDS pairs were resolved"
    elif structures["start_shift"] and not structures["stop_shift"]:
        headline = "Resolved boundary changes affect CDS starts"
    elif structures["stop_shift"] and not structures["start_shift"]:
        headline = "Resolved boundary changes affect CDS stops"
    elif structures["start_shift"] or structures["stop_shift"]:
        headline = "Both start and stop boundaries differ between annotations"
    elif changed:
        headline = "Translations differ despite identical CDS structures"
    elif limited:
        headline = "Some accepted pairs have limited protein evidence"
    else:
        headline = "Accepted CDS pairs encode identical proteins"
    descriptions = {
        "identical": "The translated residues at corresponding genomic codons are identical.",
        "terminal_only": "Only a terminal boundary differs; all shared translated residues match.",
        "initiator_recoding": "Initiator-methionine recoding explains the shared-codon mismatch; terminal length may also differ.",
        "translation_difference": "Other translation or coding-path differences; inspect genetic codes, exceptions and partial annotations.",
        "limited": "Unavailable or ambiguous translations, internal stops, or codon-mapping limits prevent residue-level interpretation.",
    }
    outcomes = ""
    for key, label in CONSEQUENCE_LABELS.items():
        count = counts.get(key, 0)
        outcomes += f'''<tr data-consequence="{key}"><td><span class="gene-dot gene-color-{key}"></span><strong>{label}</strong><p>{descriptions[key]}</p></td>
<td class="gene-result-count"><strong>{count:,}</strong><small>{percent(count, matched)}</small></td>
<td><button type="button" data-consequence="{key}" aria-controls="gene-table" {'disabled' if not count else ''}>View rows</button></td></tr>'''
    boundary_rows = "".join(f'<tr><th scope="row">{label}</th><td>{structures[key]:,}</td><td>{percent(structures[key], matched)}</td></tr>'
                            for key, label in (("exact", "Exact CDS structure"), ("start_shift", "Alternative start boundary"), ("stop_shift", "Alternative stop boundary")))
    length_summary = '<p class="note">No interpretable protein differences; length-shift summary unavailable.</p>'
    if deltas:
        length_summary = f'''<p class="gene-length-directions">B shorter: <b>{sum(delta < 0 for delta in deltas):,}</b> · B longer: <b>{sum(delta > 0 for delta in deltas):,}</b> · Same length: <b>{deltas.count(0):,}</b></p>
<p class="note">Among {len(deltas):,} interpretable protein differences, median absolute Δ: <strong>{median(abs(delta) for delta in deltas):g} aa</strong>; signed range: {min(deltas):+d} to {max(deltas):+d} aa. Δ is B minus A; limited and unresolved evidence is excluded.</p>'''
    # HTML parsers recognize </script> even inside a JSON string. Escape markup
    # before embedding; the client always renders source attributes as text.
    payload = json.dumps(gene, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    assets = files("annostat").joinpath("assets")
    css = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.css", "gene_report.css"))
    scripts = "\n".join(assets.joinpath(name).read_text(encoding="utf-8") for name in ("datatables-3.1.4.min.js", "gene_report.js"))
    license_text = escape(assets.joinpath("datatables-LICENSE.txt").read_text(encoding="utf-8"))
    return f'''<style>{css}</style><section class="gene-report" id="gene-report">
<nav class="gene-nav" aria-label="Report sections"><a href="#gene-overview">Overview</a><a href="#gene-differences">Differences</a><a href="#gene-evidence">Gene evidence</a><a href="#gene-methods">Methods</a><a href="#report-data">Download data</a></nav>
<section class="gene-section" id="gene-overview"><div class="gene-section-heading"><span>01 / Overview</span><p>Two annotations · identical genome sequences</p></div>
<div class="gene-finding"><span class="gene-eyebrow">Main finding</span><h2>{headline}</h2>
<p><strong>{matched:,} resolved one-to-one CDS pairs.</strong> Protein comparisons: {identical:,} identical, {changed:,} interpretable differences, {limited:,} limited.</p>
<p class="note">Experimental support is not assessed. These findings describe annotation differences and do not identify the biologically correct model.</p></div>
<div class="gene-counts">{cards}</div>
<div class="gene-sources">{sources}</div></section>
<section class="gene-section" id="gene-differences"><div class="gene-section-heading"><span>02 / Difference profile</span><p>Counts and percentages use accepted pairs</p></div>
<div class="gene-profile-grid"><div><h2>Translated protein differences</h2><p class="note">Denominator: {matched:,} accepted pairs. Limited evidence is retained in this denominator.</p>
<div class="table-wrap"><table class="gene-outcomes"><thead><tr><th>Finding and interpretation</th><th>Pairs / %</th><th>Evidence</th></tr></thead><tbody>{outcomes}</tbody></table></div></div>
<div><h2>CDS structure</h2><p class="note">Structural correspondence is assessed separately from protein equality.</p>
<table class="gene-boundaries"><thead><tr><th>Accepted relationship</th><th>Pairs</th><th>%</th></tr></thead><tbody>{boundary_rows}</tbody></table>
<p class="note">All-model structural agreement: <strong>{exact}</strong>. Formula: 2 × {gene["exact_model_count"]:,} shared models / {sum(gene["cds_counts"].values()):,} CDS models across both sources. Shared multiplicity is included; this denominator also contains models without an accepted counterpart.</p>
<h3 class="gene-length-heading">Protein length shifts</h3>{length_summary}
<aside class="gene-unresolved"><h3>{unresolved:,} unresolved relationship rows</h3>
<p>{gene["status_counts"].get("unmatched", 0):,} unmatched calls · {gene["status_counts"].get("candidate", 0):,} candidate rows · {gene["status_counts"].get("ambiguous", 0):,} ambiguous rows.</p>
<p class="note">Relationship rows are not distinct genes. Unmatched calls do not establish biological absence; unresolved overlaps do not establish split/merge events.</p>
<button type="button" data-consequence="unresolved" aria-controls="gene-table" {'disabled' if not unresolved else ''}>Review unresolved rows</button></aside></div></div>
<details class="gene-figure"><summary>Editable figure: protein consequence counts</summary>{plot.read_text(encoding="utf-8")}<a href="plots/gene_consequences.svg">Download SVG</a></details></section>
<section class="gene-section" id="gene-evidence"><div class="gene-section-heading"><span>03 / Gene evidence</span><p>Every relationship is available</p></div>
<h2>Inspect individual CDS models</h2><p class="note">Search a gene, CDS ID, product or source attribute. Select <strong>Inspect</strong> to connect the classification to its coordinates and translated residues. All {len(gene["relationships"]):,} relationship rows are retained.</p>
<div class="gene-controls"><label>Protein evidence <select id="gene-consequence"><option value="">All relationships</option>{''.join(f'<option value="{key}">{label}</option>' for key, label in CONSEQUENCE_LABELS.items())}<option value="unresolved">Unresolved correspondence</option></select></label>
<label>Structure <select id="gene-structure"><option value="">All structures</option>{''.join(f'<option value="{key}">{key.replace("_", " ")}</option>' for key in ("exact", "start_shift", "stop_shift", "frame_disagreement", "overlap", "unmatched"))}</select></label></div>
<noscript><p>JavaScript is required for the interactive table. The complete <a href="tables/gene_comparisons.tsv">TSV</a>, <a href="comparison.json">JSON</a>, and <a href="comparison.txt">text summary</a> remain available.</p></noscript>
<div class="table-wrap"><table id="gene-table" aria-label="Complete CDS relationships"></table></div>
<section id="gene-detail" class="gene-detail" aria-live="polite"><h3>Gene evidence</h3><p>Select a relationship to inspect its two annotation models.</p></section>
</section><section class="gene-section" id="gene-methods"><div class="gene-section-heading"><span>04 / Methods and interpretation</span><p>Descriptive measurements; no significance tests</p></div>
<h2>Read the measurements correctly</h2><p>Genomic coordinates are <strong>1-based inclusive</strong>. Segments remain in biological order; arrows show strand. Positive boundary deltas point downstream in the CDS strand direction. Protein length and terminal deltas are <strong>B minus A</strong>, in amino acids (aa).</p>
<p>{escape(gene["consequence_methods"])}</p><p>Non-gap mismatches are measured at corresponding genomic codons, with initiation-context mismatches reported separately from other mismatches. They are not DNA mutations. Gap aa counts codons present in only one translated model. Identity includes gap columns. Unavailable measurements are not zero.</p>
<p>{escape(gene["scope"])} These are descriptive annotation differences, not tests of significance, functional-impact predictions, or evidence that either annotation is correct. Percentages in the difference profile use all accepted pairs; they do not estimate annotation accuracy.</p>
<p class="downloads"><a href="comparison.txt">Interpretation summary (TXT)</a><a href="tables/gene_comparisons.tsv">All relationships and evidence (TSV)</a><a href="plots/gene_consequences.svg">Consequence counts (SVG)</a></p>
<details><summary>Bundled table-library license</summary><p>DataTables 3.1.4, locally bundled; no network access.</p><pre class="gene-license">{license_text}</pre></details>
</section></section><script id="gene-data" type="application/json">{payload}</script><script>{scripts}</script>'''


def gene_text_report(summary: dict) -> str:
    """Keep the human-readable summary short; the TSV retains every row."""
    gene = summary["gene_level_comparison"]
    matched = gene["status_counts"].get("matched", 0)
    agreement = gene["exact_agreement"]
    exact = f"{100 * agreement:.2f}%" if agreement is not None else "n/a (no CDS models)"
    labels = [profile["label"] for profile in summary["datasets"]]
    lines = [f'Annostat {summary["annostat_version"]} | CDS annotation comparison',
             f"A: {labels[0]}", f"B: {labels[1]}",
             "Scope: two annotations of identical DNA; genomic coordinates 1-based inclusive.",
             f"Resolved one-to-one CDS pairs: {matched:,}"]
    lines += [f"{label}: {gene['consequence_counts'].get(key, 0):,} / {matched:,} resolved pairs" for key, label in CONSEQUENCE_LABELS.items()]
    lines += [f"Exact structural agreement: {exact}; {gene['exact_model_count']:,} shared models (multiplicity included).",
              "Unmatched CDS calls: " + "; ".join(f"{label}: {count:,}" for label, count in gene["unmatched_cds_counts"].items()),
              f"Candidate relationship rows: {gene['status_counts'].get('candidate', 0):,}; ambiguous relationship rows: {gene['status_counts'].get('ambiguous', 0):,}.",
              "Relationship rows are not distinct genes. Unresolved pairs have no assigned protein consequence.", "",
              "Review shortlist (up to 3 examples per group; all rows remain in TSV/JSON):",
              "Within each group, measured absolute protein-length differences are largest first; unavailable lengths follow. This is not a functional-impact ranking."]
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
            for end in ("n", "c"):
                if row[f"{end}_terminal_delta_aa"]:
                    residues = row[f"{end}_terminal_residues_a"] or row[f"{end}_terminal_residues_b"]
                    lines.append(f"  {end.upper()} terminus B-A {row[f'{end}_terminal_delta_aa']:+d} aa; extra sequence in {'A' if row[f'{end}_terminal_residues_a'] else 'B'}: {residues[:60]}{'... (full sequence in TSV/JSON)' if len(residues) > 60 else ''}.")
    lines += ["", "Interpretation: differences do not establish functional effects or annotation correctness. Unmatched calls do not establish biological absence. Experimental support is not assessed by this report.",
              gene["consequence_methods"],
              "Full evidence: comparison.json; tables/gene_comparisons.tsv; interactive comparison.html; plots/gene_consequences.svg.", ""]
    return "\n".join(lines)
