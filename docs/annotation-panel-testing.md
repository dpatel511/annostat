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
cohort joins and exported SVG geometry/labels are cross-checked.
Pseudogene models without a published translation have computed-only reference
proteins and are counted separately. Pyrodigal's short-record GenBank LOCUS lines
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
