# Multi-tool public annotation audit

This optional development suite checks completed annotations on **identical DNA**.
Coordinates are **1-based inclusive**; sequence IDs are mapped only after exact
sequence/hash and length agreement. Original files and transformations are retained.

| Assembly | Organism | Annotation inputs |
| --- | --- | --- |
| GCF_000009045.1 | Bacillus subtilis | Bakta, RefSeq, Prokka, Pyrodigal |
| GCF_000829055.1 | Lactobacillus casei (source name) | Bakta, RefSeq, Prokka, Pyrodigal |
| GCF_000817005.1 | Streptococcus pneumoniae | Bakta, RefSeq, Prokka, Pyrodigal |
| GCF_000013425.1 | Staphylococcus aureus | Bakta, RefSeq, Prokka, Pyrodigal |
| GCF_000067165.1 | Sorangium cellulosum | Bakta, RefSeq, Pyrodigal |

The matrix contains 19 annotation datasets and 27 distinct tool pairs, each tested
in both directions. Sorangium has no verified Prokka input in this panel. Its
large, GC-rich assembly supplements the smaller chromosome/plasmid examples.

Bakta 1.1/database 3.0 files come from the [published benchmark examples](https://doi.org/10.5281/zenodo.5253552);
Prokka 1.14.6 outputs come from [HOMD release 11.02](https://www.homd.org/ftp/genomes/PROKKA/V11.02/).
RefSeq snapshots are identified by the exact NCBI URLs, headers and hashes.
`examples/annotation-panel-manifest.json` freezes every download and prepared input.
These historical versions test Annostat compatibility; they do not benchmark the
latest annotator releases. Source availability can change; a checksum mismatch
requires investigation, never automatic replacement of expected hashes.

[Pyrodigal](https://pyrodigal.readthedocs.io/en/stable/) 3.7.1 supplies a test-only
CDS-structure comparator, trained on each assembly using table 11. It has no
functional annotation database and shares Prodigal ancestry with the other pipelines.
Zero explicitly hypothetical labels in an input without product labels do not
establish complete functional annotation. Annotation disagreements do not identify
the correct model, prove a functional effect, or establish gene presence/absence.

From the repository in WSL/Linux, use a dedicated test environment and output directory:

```bash
python -m venv tmp/panel-venv
source tmp/panel-venv/bin/activate
python -m pip install -e . -r examples/annotation-panel-requirements.txt
python examples/prepare_annotation_panel.py --output tmp/annotation-panel
python examples/check_annotation_panel.py \
  --inputs tmp/annotation-panel --output tmp/annotation-panel/results
# Recheck existing outputs after changing the checker, without rerunning the CLI:
python examples/check_annotation_panel.py \
  --inputs tmp/annotation-panel --output tmp/annotation-panel/results --check-only
```

Preparation downloads only completed public annotations and runs the small
structural comparator. Budget roughly 2 GB for sources and the many full reports, plus the optional
browser installation; keep all data/results under ignored `tmp/`.
The runner starts at most three pairwise CLI jobs concurrently. It writes separate
stdout/stderr logs, command exit statuses, environment/source hashes and
`verification.json`. Optional dependencies belong to this test environment only.

The independent oracle imports no Annostat calculations. It joins GFF CDS models
to source GenBank/FAA records, verifies every nucleotide/protein export, enumerates
CDS overlaps independently, and checks accepted/candidate/ambiguous/unmatched
correspondence. Both directions are checked for signed deltas and codon-mapped
residue counts. GC, union coding density, feature/RNA/COG counts, codon usage,
aggregate profiles, text denominators, HTML evidence, TSV fields, group membership,
cohort joins and exported SVG geometry/labels are cross-checked. CDS row counts,
per-CDS genetic-code counts, recognized start codons and per-model outcomes are
recomputed from the source records and the oracle's own correspondence graph.
Pseudogene models without a published translation have computed-only reference
proteins and are counted separately; when a source FAA lists such a protein
(Bakta 1.12), the computed protein must equal it. Pyrodigal's short-record GenBank LOCUS lines
may produce Biopython warnings; explicit sequence, coordinate and ID joins remain
mandatory. No warning is suppressed to declare success.

The separate browser audit requires Node.js and Playwright (install under `tmp/`
to avoid adding runtime dependencies):

```bash
npm install --prefix tmp/panel-browser playwright@1.62.1
export NODE_PATH="$PWD/tmp/panel-browser/node_modules"
node tmp/panel-browser/node_modules/playwright/cli.js install chromium
node examples/check_annotation_panel_browser.cjs tmp/annotation-panel/results
```

Alternatively set `ANNOSTAT_BROWSER` to an installed Chromium/Chrome/Edge binary.
The browser checks every report at 1440 and 430 pixels, category filters,
representative evidence types, neighbouring context, SVG downloads, displayed
counts, offline operation and clipping. An explicit `en-IN` browser locale
checks that coordinate/count grouping still matches the exported reports. It also checks every standalone plot for
text outside the exported canvas. Screenshots and results are retained under
`results/browser/`. Browser checks supplement scientific checks; a rendered SVG
alone does not establish publication quality or biological correctness.

The 2026-10-31 local run passed all 99 CLI checks, 155 unit tests and 314 browser/
SVG checks (79 reports at two widths and 156 plots). It cross-checked 76,407
source protein records; 641 pseudogene models used computed-only references.
All real-data inputs use table 11; other genetic codes and fragmented draft
assemblies are outside this panel. Detailed evidence stays in the chosen output
directory rather than being distributed with the package.

## Description-summary audit

The panel oracle also checks source-derived product-status assignments,
identical-protein eligibility, all nine matrix cells, gene-symbol missingness,
distinct unresolved models, TSV fields and text denominators. This uses
independently parsed source attributes and source protein/codon correspondence,
rather than importing Annostat's summary calculations.

Run the small constructed cases without downloading public inputs:

```bash
python examples/check_description_summary.py --output tmp/description-cases
node examples/check_annotation_panel_browser.cjs tmp/description-cases
```

The 29 examples run in both directions and cover both strands, nonzero initial
phases, multipart/circular CDSs, tables 4/11/25, partial starts, supported and
unapplied translation exceptions, unknown/internal-stop residues, changed
boundaries, inherited metadata, missing/uncertain labels, escaped Unicode/HTML,
duplicate structures, ambiguous overlaps, frame candidates and empty sources.
Their expected outcomes are declared fixtures, not inferred from the production
description-summary helper. The browser step uses the optional setup above.

Use fresh output directories after a schema/feature change. `--check-only`
checks existing outputs against the current contract; it does not upgrade old
reports.

On 2026-11-13 the extended local audit checked the same frozen five-assembly
panel: **99 real-data CLI commands**, **19 analyses**, **27 pairs in both
directions**, **76,407 published source proteins** and **641 computed-only
pseudogene models**. The 29 constructed examples also passed on WSL and native
Windows in both directions. Four public-input filtering/profiling controls
verified selected IDs/sequences and unchanged full-analysis metrics, including
an empty COG-required export from a source without category fields.

All **171 unit/subprocess tests** passed on Python 3.10.22 (Biopython 1.88,
Matplotlib 3.10.9) and Windows Python 3.12.14 (Biopython 1.88, Matplotlib 3.11.2,
UTF-8 mode disabled, CP1252 locale). The audit reproduced and fixed a terminal
summary crash on nonrepresentable labels/paths. Its regression and installed-wheel
checks force strict CP1252 output while verifying original Unicode in UTF-8
reports. Nine deliberately corrupted result/export controls were rejected,
and ten invalid-input/validation-exit-policy checks passed.

Offline rendered checks passed for **199 HTML reports at two widths and 399
standalone SVGs**: 797 browser/figure cases across real, synthetic, Windows and
filtered outputs. Representative desktop/mobile screenshots were inspected.
Detailed locus/source tables may scroll horizontally on narrow screens;
no page-level overflow or exported SVG text clipping was detected.

The final wheel was checked outside the checkout on Python 3.10 and native
Windows; all 26 runtime files matched the wheel and checkout. No annotation
database or new genome data was downloaded. These are computational/source
consistency checks, not experimental biological validation or a ranking of
current annotation-tool releases. Exact logs, source hashes, original failure
and final wheel checks remain under the ignored audit directory.

## Current-version audit (2026-11)

`examples/annotation-audit-2026-11-manifest.json` pins a second panel of eight
NCBI reference assemblies. The five panel assemblies have three annotations each:
the RefSeq records NCBI served on 2026-11-21, the Bakta 1.1 examples above and
Bakta 1.12.2 (database 6.0, full). Three assemblies added on 2026-11-22 cover NCBI
table 4 and draft assemblies with RefSeq and Bakta 1.12.2. It does not replace the
frozen panel. Each RefSeq version is read from its record:

| Assembly | Coverage | RefSeq annotation |
| --- | --- | --- |
| GCF_000009045.1 | complete, table 11 | no PGAP annotation comment |
| GCF_000829055.1 | complete, table 11 | PGAP 6.10, 2026-03-30 |
| GCF_000817005.1 | complete, table 11 | PGAP 6.3, 2022-11-11 |
| GCF_000013425.1 | complete, table 11 | no PGAP annotation comment |
| GCF_000067165.1 | complete, table 11 | PGAP 6.10, 2026-03-31; 109 product and 38 inference values changed since the frozen snapshot with no coordinate or date change |
| GCF_040556925.1 (*Mycoplasmoides genitalium* G37) | complete, table 4 | PGAP 6.9, 2025-03-31 |
| GCF_025234795.1 (*Exiguobacterium mexicanum*) | draft, 135 contigs | PGAP 6.11, 2026-08-03 |
| GCF_006539825.1 (*Corynebacterium variabile*) | draft, 68 contigs | PGAP 6.11, 2026-09-10 |

Bakta 1.12.2 was run through the [Bakta Web](https://bakta.computational.bio) API
(Bakta 1.12.2, database 6.0, backend 1.1.0) on each RefSeq genomic FASTA with
original contig headers, the assembly's table and the genus/species/strain recorded
in the manifest; only complete assemblies were submitted as complete genomes. Each
result records its Bakta and database versions in its JSON and GFF3 header. Bakta
1.12.2 GFF3 declares no `transl_table` (its GenBank file does), so the manifest
records table 4 for that annotation and the checker passes it as `--genetic-code`.
Without it, Annostat translates the *M. genitalium* Bakta models with table 11:
372 of 522 proteins contain internal stops, reported as QC warnings, and the
checker rejects the exports. Bakta Web deletes jobs after 30 days and the outputs
are not deposited, so the manifest pins them by hash; a new run may not reproduce
identical bytes and then needs a reviewed new manifest. Supply the outputs as
`<accession>/bakta.{fna,gff3,gbff,faa}`:

```bash
python examples/prepare_annotation_panel.py --manifest examples/annotation-audit-2026-11-manifest.json \
  --supplied path/to/bakta-outputs --output tmp/current-audit
python examples/check_annotation_panel.py --manifest examples/annotation-audit-2026-11-manifest.json \
  --inputs tmp/current-audit --output tmp/current-audit/results
```

NCBI translates an unambiguous incomplete final codon of a 3'-partial CDS; Annostat
excludes incomplete trailing bases. The oracle therefore compares such models
without that one residue (18 RefSeq CDSs in the two drafts, counted per source as
`partial_final_codon_residues`) after checking the GenBank translation against the
FAA unchanged.

The 2026-11-22 run with Annostat 1.3.0 plus the circular-outcome fix (WSL, Python
3.10) passed **21 analyses** and **18 pairs in both directions**, cross-checking
**77,336 published source proteins** with 738 computed-only pseudogene models.
Bakta 1.1 and 1.12.2 have identical structures and proteins for every accepted
pair. Seven corrupted outputs, one per new oracle check, and the table-4 case
without its recorded table were rejected. The offline browser audit passed 277
report and figure checks.

Limits: table 4 rests on one small genome and drafts on two PGAP 6.11 references;
table 25 and fungal annotations are not covered. No public pair contains an
identical duplicate structure or a same-ends internal change, so those outcomes
rest on unit tests. A constructed circular case showed that 1.3.0 classed a
same-ends internal change on a circular sequence as a partial overlap; this is
fixed after 1.3.0.

## Spliced fungal (table 1) checks

`examples/check_fungal_comparison.py` independently checks `analyze` and
`compare --gene-level` outputs for gene → mRNA → CDS annotations. It imports no
Annostat calculations. It recounts genes, mRNAs, CDS models and CDS rows,
retranslates every model with Biopython `CompoundLocation` extraction, and
optionally joins source proteins by `protein_id`, transcript `Parent` or CDS ID.
For comparisons it recomputes exact structural matches and requires accepted
pairs with internal stops to be `limited`. With `--expect-self` it also requires
every unpaired model to belong to a CDS-identical isoform group.

```bash
python examples/check_fungal_comparison.py --genetic-code 1 \
  --analysis results/a genome.fa a.gff3 --source-proteins a.proteins.fa \
  --comparison results/self A genome.fa a.gff3 B genome.fa a.gff3 --expect-self
```

Documented export conventions are classified rather than treated as failures:
GenBank translations of an incomplete final codon, source resolution of codons
with non-ACGT bases (Annostat keeps `X`), and `transl_except` recoding, which
is checked against source proteins only.
