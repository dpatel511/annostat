# Annostat

[![PyPI](https://img.shields.io/pypi/v/annostat)](https://pypi.org/project/annostat/)
[![Python](https://img.shields.io/pypi/pyversions/annostat)](https://pypi.org/project/annostat/)
[![CI](https://github.com/dpatel511/annostat/actions/workflows/ci.yml/badge.svg)](https://github.com/dpatel511/annostat/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dpatel511/annostat/blob/main/LICENSE)

**Annostat** is a command-line toolkit for validating,
analyzing, comparing, and summarizing bacterial genome annotations from matching
GFF3 and FASTA files. It produces reproducible machine-readable results,
self-contained HTML reports, CDS sequence exports, and editable SVG figures.

> Annostat evaluates annotation structure and descriptive quality indicators. It
> does not replace an annotation pipeline, infer orthology, calculate ANI, or
> establish biological correctness from sequence alone.

This README describes **Annostat 1.1.0**, which adds opt-in CDS-by-CDS comparison
of two annotations on the same unchanged genome. Check `annostat --version`:
package indexes may provide an earlier release, and unreleased changes require
a source checkout containing those changes.

## Contents

- [Capabilities](#capabilities)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command reference](#command-reference)
- [Common workflows](#common-workflows)
- [Inputs and scientific scope](#inputs-and-scientific-scope)
- [Outputs](#outputs)
- [Reproducibility and interpretation](#reproducibility-and-interpretation)
- [Help and troubleshooting](#help-and-troubleshooting)
- [Testing and development](#testing-and-development)
- [Scientific basis](#scientific-basis)

## Capabilities

- Counts CDS, RNA types, all feature types, hypothetical proteins, gene annotations,
  and CDS with COG categories
- Exports a compact feature overview as CSV or TSV
- Extracts every CDS on the correct strand as nucleotide and amino-acid multi-FASTA
- Calculates codon usage percentages, start-codon counts, and individual COG category counts
- Reports GC content, coding density, CDS overlaps and containment, pseudogenes,
  partial features, structural-RNA completeness, frame issues, and ambiguous bases
- Supports optional CDS exports filtered by length, COG annotation, and
  hypothetical-protein status without changing the complete analysis
- Builds a self-contained offline HTML report with provenance and performance timings
- Saves three focused, editable SVG charts for COG categories, CDS lengths,
  and table-aware start codons; the COG chart is omitted when the source has no
  COG data
- Streams CDS processing and FASTA generation to reduce runtime and peak memory
- Compares two or more local or NCBI-hosted annotations using normalized genome,
  annotation-completeness, COG, start-codon, and QC profiles
- Reports pairwise Jensen-Shannon profile distances and exact gene-symbol/COG-ID
  overlap without presenting annotation labels as inferred orthologs
- Optionally compares individual CDS models from two annotations of identical
  genome sequences with `compare --gene-level`
- Downloads versioned RefSeq or GenBank assemblies through the official NCBI
  Datasets command-line tool when network connectivity is requested
- Separates deterministic FASTA/GFF3 validation from biological inspection, with
  stable rule IDs, severity levels, schema/ruleset versions, and source citations
- Records SHA-256 hashes for reproducibility and writes canonical JSON plus TSV
  validation results suitable for CI pipelines
- Correctly joins repeated-ID multipart CDS features in biological order on both
  strands and validates phase continuity before translation
- Supports NCBI prokaryotic translation tables 4, 11, and 25, automatically
  selecting a consistent CDS or `region` `transl_table` declaration when present
- Applies position-specific NCBI `transl_except` annotations for selenocysteine,
  pyrrolysine, explicitly annotated termination codons, and `OTHER` residues
  represented conservatively as `X`
- Aggregates completed inspections into a script-free cohort HTML report,
  normalized TSV matrix, and machine-readable JSON

### Requirements

| Component | Status | Purpose |
|---|---|---|
| Python 3.10 or newer | Required | Runs Annostat and its test suite |
| Matplotlib 3.10 or newer | Required | Generates headless, editable SVG plots |
| Biopython 1.85 or newer | Required | Provides vetted translation-table metadata and sequence-alignment utilities |
| Genome FASTA (`.fna`, `.fa`, or `.fasta`) | Required | Supplies the genome sequences |
| Matching GFF3 annotation (`.gff` or `.gff3`) | Required | Supplies feature coordinates and annotations |
| Matching FASTA/GFF3 sequence identifiers | Required | Connects every annotated feature to its sequence |
| NCBI Datasets CLI | Optional | Downloads GCF/GCA assemblies for `fetch` or `compare --reference` |
| Internet connection | Optional | Needed only when downloading NCBI data |

### Core outputs

Annostat provides the following core analysis outputs:

- Counts CDS, RNA features by type, hypothetical CDS, CDS with gene names, CDS
  with COG categories, and other feature types.
- Writes a CSV or TSV feature overview containing ID, type, start, stop, gene,
  product, and strand, with missing values left empty.
- Exports all CDS as nucleotide and translated amino-acid multi-FASTA files.
- Calculates pooled codon usage percentages, start-codon counts, and detailed
  COG-category statistics, including multi-category assignments.
- Generates at least two meaningful scientific plots; a normal run produces CDS
  length, start-codon, and—when available—COG-category SVG figures.
- Provides an `argparse` command-line interface and modular Python structure with
  documented functions, classes, and methods.

Filtering, quality control, HTML reports, performance profiling, comparative
analysis, cohort summaries, and NCBI connectivity extend the core workflow.

## Inputs and scientific scope

- GFF3 records must contain nine tab-separated fields with 1-based inclusive
  coordinates. The FASTA identifier is the first whitespace-delimited header token.
- CDS strand and phase are honored. Circular-origin CDS coordinates extending past
  the sequence length are supported when a `region` feature declares
  `Is_circular=true`.
- Repeated CDS IDs are interpreted as one multipart biological feature. Segments
  must use the same sequence ID, strand, and Parent; phase continuity is checked,
  and segments are joined in biological order before translation.
- COG coverage describes categories actually present in the annotation source. An
  annotation with no COG category fields is reported as unavailable, not as 0%.
- Annostat uses the consistent CDS or `region` `transl_table` declared in GFF3.
  When the annotation declares no table, it falls back to bacterial table 11. `--genetic-code` is an
  explicit override for annotations without that metadata; an override that
  conflicts with a declared table is rejected instead of silently mistranslating
  CDS features. Annostat never guesses a genetic code from sequence composition.
- Genome GC percentage uses canonical A, C, G, and T bases; ambiguous bases are
  excluded from its denominator. Assembly length and coding density still use
  the complete sequence length.

## Installation

Annostat requires Python 3.10 or newer. Installation automatically includes
Matplotlib 3.10 or newer for headless SVG plot generation.

### Install from PyPI

Create an isolated environment and install the current release:

Linux or WSL:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install annostat
```

Windows PowerShell:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install annostat
```

Verify the installation:

```bash
annostat --version
annostat --help
```

Expected output for a 1.1.0 installation:

```text
Annostat 1.1.0
```

### Install from source

Use a source checkout for development or an unreleased fix:

```bash
git clone https://github.com/dpatel511/annostat.git
cd annostat
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

On Windows, activate with `.\.venv\Scripts\Activate.ps1`. If the installed
executable is unavailable, the equivalent module form is `python -m annostat`.

### Bioconda

Annostat is available through [Bioconda](https://bioconda.github.io/recipes/annostat/README.html).
On Linux, WSL, or macOS, create an environment using the
[recommended channel order and strict priority](https://bioconda.github.io/#with-conda):

```bash
conda create -n annostat --channel conda-forge --channel bioconda \
  --strict-channel-priority annostat
conda activate annostat
annostat --version
```

Bioconda releases may lag PyPI or GitHub. The `--gene-level` option requires
Annostat 1.1.0 or newer. For native Windows, use the PyPI instructions above.

## Quick start

### 1. Check the inputs

Confirm that you have one genome FASTA and its matching GFF3 annotation. Sequence
IDs must agree—for example, a FASTA record named `chromosome` must be referenced
as `chromosome` in the first GFF3 column.

### 2. Validate the paired files

Run deterministic integrity checks before biological interpretation:

```bash
annostat validate \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o validation
```

Example output for the included dataset:

```text
Annostat 1.1.0 | annotation validation PASS
  No validation issues found; no output files were written.
```

When no issues are found, the standalone command does not create an output
directory or files. If it finds any warning or error, it writes
`validation.json` and `validation.tsv`. Exit status `0` means the configured
failure threshold was not reached; status `1` means it was reached. Use
`--fail-on warning` for strict CI, or `--fail-on never` when the report should
never stop a pipeline. A validation pass means the paired files are internally
interpretable; it does not prove that every biological annotation is correct.

### 3. Analyze one annotation

```bash
annostat analyze \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

Here `-f` selects the required FASTA, `-g` selects the required GFF3, and `-o`
chooses the output directory. Add `--table-format tsv` when a TSV overview is
preferred over CSV. Translation-table selection is automatic when the
GFF3 CDS or `region` rows contain `transl_table`; otherwise table 11 is used.
The included example produces a 5,076,188 bp analysis with 4,295 CDS. For a local
annotation without translation-table metadata, select a supported table explicitly
with `--genetic-code 4`, `--genetic-code 11`, or `--genetic-code 25`.

### 4. Review the result package

Open `results/report.html` in a browser first. Then use:

- `summary.json` for the complete machine-readable summary;
- `tables/` for feature, codon, COG, start-codon, and QC data;
- `sequences/` for nucleotide and protein CDS FASTA files; and
- `plots/` for editable SVG figures.

### 5. Summarize a cohort

After inspecting several genomes into separate result directories, list those
directories explicitly or pass their common parent directory:

```bash
annostat summarize results/sample-a results/sample-b -o cohort
```

The output contains `cohort.html`, `cohort.tsv`, and `cohort.json`. Samples are
ordered deterministically. Missing source-dependent values such as COG coverage
remain `NA`/`null`, never biological zero. Multiple source Annostat versions are
retained and flagged in the report.

## Command reference

| Command | Purpose | Primary output |
|---|---|---|
| `annostat analyze` | Analyze one matching FASTA/GFF3 pair | HTML, JSON, tables, CDS FASTA, SVG |
| `annostat validate` | Check deterministic structural and cross-file integrity | JSON and TSV validation findings |
| `annostat summarize` | Aggregate completed Annostat analyses | Cohort HTML, TSV, and JSON |
| `annostat compare` | Compare two or more local or NCBI annotations | Comparative HTML, JSON, TSV, and SVG |
| `annostat compare --gene-level` | Compare CDS models in exactly two annotations of the same unchanged genome | Adds CDS relationships to the comparison HTML, JSON, and TSV |
| `annostat fetch` | Download versioned NCBI assemblies | FASTA, GFF3, and assembly metadata |

Run `annostat COMMAND --help` for required arguments, defaults, examples, and
output behavior. `annostat inspect` is a fully supported alias for the same
single-annotation workflow. The older `annostat -f ... -g ...` form remains
compatible in the 1.x series but is deprecated; new scripts should use an
explicit command.

## Common workflows

The following sections show filtered exports, comparative analysis, NCBI
references, performance profiling, and complete output descriptions.

Without installing, module execution remains available from the repository root:
`python -m annostat --help`.

Use `--table-format tsv` for a tab-separated feature overview. Add `--profile`
to measure peak Python memory and print per-stage timings.

### Filtered CDS exports

Filters create an additional `filtered/` directory while retaining every normal
analysis file. Criteria are combined, so a CDS must satisfy all enabled options:

```bash
annostat analyze \
  -f genome.fna \
  -g annotations.gff3 \
  -o results \
  --min-cds-length 300 \
  --require-cog \
  --exclude-hypothetical
```

Available criteria are `--min-cds-length`, `--max-cds-length`, `--require-cog`,
and `--exclude-hypothetical`.

### Comparative analysis

Compare two or more labelled local genomes with repeatable `--genome` options:

```bash
annostat compare \
  --genome genome-a data/GCF_000007145.1.fna data/GCF_000007145.1.gff3 \
  --genome genome-b data/GCF_001050915.2.fna data/GCF_001050915.2.gff3 \
  --output comparison
```

Each genome is structurally validated and independently uses its declared
`transl_table` or the table-11 fallback. For an annotation without
`transl_table`, a per-genome override uses
the same label supplied to `--genome` (or the reference accession):

```bash
annostat compare \
  --genome sr1 sr1.fna sr1.gff3 \
  --genome ecoli ecoli.fna ecoli.gff3 \
  --genetic-code sr1 25 \
  --output comparison
```

The comparison package is intentionally compact:

```text
comparison/
|-- comparison.html
|-- comparison.json
|-- tables/
|   |-- dataset_metrics.tsv
|   `-- pairwise_comparisons.tsv
`-- plots/
    |-- dataset_overview.svg
    `-- cog_comparison.svg (only when every input supplies COG categories)
```

`dataset_metrics.tsv` provides raw genome/feature counts together with GC,
coding density, CDS per Mb, annotation coverage, start-codon usage, and QC
warnings. `pairwise_comparisons.tsv` contains metric deltas, percentage deltas,
Jensen-Shannon distances for COG and start-codon profiles, and Jaccard overlap
for exact gene symbols and explicit COG identifiers.

For two datasets, the COG figure shows percentage-point differences. For three
or more datasets, it becomes a dataset-by-category heatmap. Exact gene-symbol
overlap describes annotation labels only; it is not orthology, ANI, phylogeny,
or a core/accessory pangenome analysis.

### Compare individual CDS annotations on the same genome

Use the existing comparison command with `--gene-level` to review changes
between annotation pipelines or annotation versions on the same genome:

```bash
annostat compare --gene-level \
  --genome A genome.fna annotation-a.gff3 \
  --genome B genome.fna annotation-b.gff3 \
  -o comparison
```

This mode requires exactly two annotations, identical FASTA sequence IDs and
bases, and explicit `+`/`-` CDS strands. FASTA descriptions, wrapping, record
order, and letter case may differ. Both GFF3 files must use the same unchanged
assembly coordinates; supplying a FASTA does not prove a GFF3 file's original
assembly. Rotated circular genomes and changed assemblies require coordinate
reconciliation before comparison.

The normal outputs remain available. `tables/gene_comparisons.tsv`,
`comparison.json`, and a new section in `comparison.html` additionally report
exact CDS models, unique start/stop alternatives, frame disagreements, other
overlap candidates, ambiguous relationships, unmatched calls, annotation
attribute changes, and translated-protein consequences. CDS sequences are
processed once per annotation; coordinates, metadata, hashes, and translated
proteins are retained for matching and consequence calculations, not CDS DNA.

All reported coordinates are **1-based inclusive**. Multipart CDSs count once,
and circular-origin blocks retain their biological order. Positive boundary
deltas point downstream in the CDS strand direction; circular deltas are only
reported for exact matches or supported boundary alternatives. Names and locus
tags do not establish correspondence. A counterpart must overlap on the same
sequence and strand. Exact models anchor matching, followed by unique,
frame-compatible boundary alternatives. Other overlaps remain candidates or
ambiguous relationships; many-to-one relationships are not automatically
declared split/merge events. RNA features and gene-only pseudogene records are
outside the CDS matching scope; unambiguous parent-gene attributes are retained.

Exact structural agreement is `2 * E / (N_A + N_B)`, where `N_A` and `N_B` count
biological CDS models and `E` sums the shared multiplicity of each exact
structure (ordered coding blocks and phase). Duplicate models can therefore
agree structurally while their ID correspondence remains ambiguous. Empty CDS
sets have unavailable agreement, not 100%. Status/relationship counts count
table rows, not distinct genes. Unmatched calls do not demonstrate biological
absence, and changed descriptions do not establish conflicting functions.
The HTML contains every relationship in a searchable, sortable, paginated table,
with protein-evidence and structural filters. Selecting a row opens the two
source models, an editable/downloadable locus SVG, changed residues, and a
protein-alignment window centered on the first difference. Full alignments are
expandable. DataTables 3.1.4 is bundled under its MIT license; the self-contained
HTML makes no CDN requests. JavaScript is needed for the table; TSV, JSON, and
text outputs remain usable without it.

`comparison.txt` provides a concise interpretation summary and up to three
examples per review group: terminal boundary differences, translation context,
limited protein evidence, and unresolved correspondence. Measured absolute length
changes are largest first within each group; unavailable lengths follow. This is
not a functional-impact ranking. The HTML opens with the main findings, per-source
pairing coverage, and separate protein/structural breakdowns. Percentages name
their denominators; all-model structural agreement includes its formula.
Protein-length direction, median absolute change, and signed range summarize
interpretable differences only, excluding limited and unresolved evidence.
Each protein category links directly to its evidence rows. Gene details separate
the interpretation, boundary measurements, locus diagram, and alignment, with
full source context expandable. Experimental support is not assessed.
`plots/gene_consequences.svg` counts protein consequences among resolved pairs
only. Re-running without `--gene-level` removes these two generated files and
the gene comparison table, and returns to the aggregate report.

Protein evidence has five categories: `identical`, `terminal_only`,
`initiator_recoding`, `translation_difference`, and `limited`. Unaccepted
correspondence is `unresolved` and excluded from protein-consequence counts.
Accepted pairs use **ordered genomic codon coordinates**, preserving initial
phase, strand, multipart order, and circular-origin wrapping. No protein-score
optimization is used to move residues away from known genomic correspondence.
An alternative start codon can translate as M at initiation but as V/L/I within
the longer CDS; this is reported explicitly as initiation context, not an
internal DNA mutation. Different tables/exceptions can change the predicted
translation even when the CDS coordinates are identical.

Identity is identical residues divided by **all alignment columns, including
gaps**; substitutions count non-gap mismatches, and `gap_aa` counts columns
containing a gap. Length/terminal deltas are **B minus A** in amino acids.
Empty translations, unknown/internal-stop residues, and missing, repeated, or
unordered codon correspondence receive limited evidence, not fabricated zero
differences. Known codon mappings do not need quadratic dynamic programming or
its resource cutoff. Terminal-only classification requires a compatible changed
CDS boundary: terminal stop reassignment at unchanged coordinates is a translation
context difference. These calculations do not predict functional effects, select
a correct annotation, or provide statistical tests.

#### Use the report for annotation review

For a runnable WSL demonstration, verified real-data examples, and an explanation
of what to show, see the [gene-level walkthrough](docs/gene-level-demo.md).

Start with `comparison.html` for resolved protein consequences, then filter and
inspect individual relationships. Aggregate profiles are secondary. Use
`comparison.txt` for a short summary and `tables/gene_comparisons.tsv` for all
side-by-side evidence:

| Finding | Useful follow-up | Interpretation limit |
|---|---|---|
| Exact structure and identical translated protein | Identify unchanged CDS models | Agreement is not independent evidence of correctness |
| `start_shift` or `stop_shift` | Inspect the boundary and protein-length difference | Neither annotation is automatically preferred |
| `frame_disagreement` | Inspect phase, coordinates, and source exception annotations | Not proof of a biological frameshift |
| Identical protein with a changed product description | Separate label changes from protein-sequence changes | Different descriptions do not prove different functions |
| `unmatched`, `candidate`, or `ambiguous` correspondence | Reconcile calls using the source annotations and external evidence | Not proof of gene absence or a split/merge event |

The complete table includes:

| Column(s) | Meaning |
|---|---|
| `dataset_a`, `dataset_b`, `id_a`, `id_b` | Source labels and CDS identifiers |
| `seqid`, `strand`, `segments_a`, `segments_b`, `phase_a`, `phase_b` | Location and coding structure; segment coordinates are 1-based inclusive |
| `cds_length_a`, `cds_length_b`, `protein_length_a`, `protein_length_b` | CDS length in nucleotides and translated-protein length in amino acids |
| `gene_a`, `gene_b`, `product_a`, `product_b`, `dbxref_a`, `dbxref_b` | Side-by-side names, descriptions, and database cross-references when supplied |
| `overlap_bp`, `start_delta_bp`, `stop_delta_bp` | Shared physical bases and strand-aware boundary shifts from A to B |
| `protein_identical`, `changed_attributes` | Equality of Annostat-translated protein sequences and the names of differing annotation attributes |
| `start_codon_a`, `start_codon_b`, `genetic_code_a`, `genetic_code_b` | First translated codon and translation table; not proof of initiation |
| `protein_consequence`, `consequence_note` | Descriptive protein category and reason; unresolved correspondence is not a protein comparison |
| `protein_length_delta_aa`, `n_terminal_delta_aa`, `c_terminal_delta_aa` | B minus A in aa; terminal deltas only when measurable, otherwise null |
| `identical_aa`, `substitution_count`, `gap_aa`, `alignment_identity_percent` | Explicit aligned counts and identity including gap columns; null if not assigned |
| `shared_codon_count`, `initiator_mismatch_count`, `internal_mismatch_count` | Shared genomic codons, initiation-specific mismatches, and other codon-mapped mismatches; not DNA mutations |
| `n_terminal_residues_a/b`, `c_terminal_residues_a/b` | Extra terminal residues in the indicated predicted protein; their cause is distinguished by the consequence category |
| `alignment_method`, `alignment_ambiguous`, `alignment_a`, `alignment_b` | Method and aligned sequences for changed, resolvable pairs; unchanged sequences are not duplicated |

Gene-level JSON schema **1.2** adds these fields, `consequence_counts`,
`unmatched_cds_counts`, and method details. Existing matching/status rules and
columns remain unchanged. In TSV, null values are blank; in JSON they are `null`.
The legacy `protein_identical` field records literal translation equality, even
for strings containing unknown residues; `protein_consequence=limited` prevents
that equality from being presented as reliable residue-level evidence.

Analysis, comparison, and cohort JSON schemas are **1.1**. The
`calculation_policy` marker records the shared normalization/metric policy.
Unambiguous parent-gene attributes are applied consistently before translation,
QC, and counting; explicit CDS attributes take precedence and original input
files are not edited. Feature tables are normalized views, not literal GFF dumps.
Blank gene values are not counted as names. COG-ID overlap is evaluated
independently of letter-category availability. GC with no A/C/G/T bases and
CDS-based percentages/medians with no CDSs are unavailable, not zero. Comparison
TSVs use blank nulls; cohort tables retain their `NA` representation. Cohorts
disclose mixed calculation policies, including legacy inputs; regenerate old
analyses consistently before interpreting cross-policy differences.

Genetic-code selections, pseudogene/partial markers, and translation-exception
attributes are also retained. Missing values are empty in TSV and `null` in JSON;
an unmatched row has no protein-equality result, not a value of `false`.
The JSON data is under `gene_level_comparison`, with CDS counts, exact-model
agreement, status/relationship row counts, and the full relationship list.

The HTML preview does not yet display all these columns or rank findings by
biological importance. Use the TSV or JSON for complete review. Neither the
agreement percentage nor the relationship labels are accuracy estimates,
statistical significance tests, or evidence of gene gain/loss between organisms.

### Compare against external NCBI data

NCBI access is optional and uses NCBI's maintained `datasets` executable. Follow
the official [installation instructions](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/),
then add an external reference directly to a local comparison:

```bash
annostat compare \
  --genome local-sample sample.fna sample.gff3 \
  --reference GCF_000007145.1 \
  --output comparison
```

`--reference` is repeatable, so a comparison may contain only external
assemblies or any mixture of local and NCBI inputs. Annostat downloads genome
FASTA, GFF3, and assembly metadata under `comparison/external_inputs/`, then
runs every input through the same comparison engine:

```bash
annostat compare \
  --reference GCF_000007145.1 \
  --reference GCF_001050915.2 \
  --output comparison
```

The report identifies organisms and whether each pair belongs to the same or
different species when the input metadata makes that determination possible.
It never assumes that labelled datasets are strains. NCBI assembly names,
taxonomy identifiers, annotation providers, and release metadata are retained
in `comparison.json`. Versioned GCF/GCA accessions in the HTML provenance table
link to their NCBI Datasets assembly pages.

COG results are source-dependent. Bakta annotations commonly include COG
categories, while current NCBI RefSeq GFF3 packages may not. When a source does
not supply COG fields, Annostat reports COG coverage and pairwise COG statistics
as `Not available` and omits the COG comparison plot; missing data is never
presented as biological zero coverage.

For users who only need the source files, the standalone downloader remains
available:

```bash
annostat fetch GCF_000007145.1 GCF_001050915.2 --output ncbi_data
```

The `NCBI_API_KEY` environment variable is honored by the official client. Local
analysis never requires a network connection, and Annostat does not transmit
local FASTA or GFF3 files.

## Outputs

### Single-annotation analysis

The `annostat analyze` output directory contains:

```text
results/
├── report.html
├── summary.json
├── tables/
│   ├── features.csv (or features.tsv)
│   ├── codon_usage.csv
│   ├── start_codons.csv
│   ├── cog_categories.csv
│   └── annotation_issues.csv
├── validation/
│   ├── validation.json
│   └── validation.tsv (when findings are present)
├── sequences/
│   ├── cds_nucleotide.fasta
│   └── cds_protein.fasta
└── plots/
    ├── cog_categories.svg (when COG categories are available)
    ├── cds_lengths.svg
    └── start_codons.svg
```

When filters are active, `filtered/` additionally contains the selected feature
table and paired nucleotide/protein FASTA files.

`report.html` works without a server or internet connection and embeds the SVG
figures directly. The separate SVG files remain available for editing and use in
documents. Codon usage remains available as the required percentage table; a
dense 64-codon heatmap is intentionally not generated by default.

`annotation_issues.csv` is a review aid rather than an automatic pass/fail
verdict. Short CDS overlaps are common in compact bacterial genomes, so ordinary
overlaps and adjacent duplicate labels are informational. Containment, incomplete
frames, missing structural RNA classes, internal stops, and ambiguous bases are
reported as warnings that merit inspection in their biological context.

### Comparison and cohort reports

`annostat compare` writes `comparison.html`, `comparison.json`,
`tables/dataset_metrics.tsv`, `tables/pairwise_comparisons.tsv`, and SVG plots
as described in [Comparative analysis](#comparative-analysis). With
`--gene-level`, it additionally writes `tables/gene_comparisons.tsv`,
`comparison.txt`, and `plots/gene_consequences.svg`, and adds
the interactive CDS evidence section and JSON results described
[above](#compare-individual-cds-annotations-on-the-same-genome). The gene-level
option does not export a corrected or consensus annotation.

`annostat summarize` writes `cohort.html`, `cohort.tsv`, and `cohort.json` from
completed single-annotation analyses. It does not aggregate gene-level
relationship tables.

## Reproducibility and interpretation

### Validation versus analysis

`annostat validate` answers whether the files can be interpreted consistently:
GFF3 syntax, FASTA symbols, matching sequence IDs, coordinate bounds, Parent/ID
relationships, multipart-feature consistency, CDS phase continuity, and valid,
consistent CDS translation-table declarations. These checks are deterministic
and may safely control a CI exit status.

`annostat analyze` first requires structural validation to pass, then reports
biological review candidates such as internal stops, missing terminal codons,
ambiguous CDS bases, overlaps, partial features, pseudogenes, duplicated adjacent
labels, and structural-RNA coverage. These findings are context-dependent and do
not automatically make an annotation invalid. Boundary-specific partial markers
suppress only the affected start- or stop-codon assumption. Explicit `pseudo`,
`transl_except`, and `exception` annotations suppress checks whose assumptions do
not apply.

Every validation result embeds its rule definitions, rationale, source URL,
ruleset version, schema version, Annostat version, and input SHA-256 hashes. JSON
and TSV findings use canonical ordering and contain no timestamp, so repeated runs
on identical paths and bytes produce byte-identical validation artifacts. A
`scientific_fingerprint` in both validation and analysis JSON excludes file
paths and performance timings, allowing scientific results to be compared across
machines and output directories.

## Help and troubleshooting

Start with the command-specific help rather than guessing an option:

```bash
annostat --help
annostat analyze --help
annostat validate --help
annostat summarize --help
annostat compare --help
annostat fetch --help
```

Common problems:

- **`annostat` is not found:** activate the environment where it was installed,
  or run `python -m annostat --help`.
- **The wrong Python starts:** use `python3 -m pip` on Linux/WSL or `py -3 -m pip`
  on Windows, then invoke Annostat from the same environment.
- **FASTA and GFF3 identifiers do not match:** compare the FASTA header token with
  the first GFF3 column; Annostat intentionally rejects ambiguous pairings.
- **COG values are unavailable:** the annotation source did not provide COG
  categories. Annostat does not reinterpret missing source data as zero.
- **`--gene-level` is unrecognized:** check `annostat --version`; this option
  requires 1.1.0 or newer and belongs to `annostat compare`, not `analyze`.
- **Gene-level comparison rejects the genomes:** it requires exactly two
  annotations with identical FASTA identifiers and bases on unchanged assembly
  coordinates. Different strains, changed assemblies, or rotated circular
  sequences are not supported without prior coordinate reconciliation. Do not
  rename identifiers merely to bypass this safeguard.
- **NCBI download fails:** confirm that the official NCBI Datasets CLI is installed
  and available on `PATH`; local analysis itself never needs network access.

For reproducible bug reports, include `annostat --version`, the complete command,
the terminal error, and a minimal non-sensitive input pair when possible. Report
problems through the [GitHub issue tracker](https://github.com/dpatel511/annostat/issues).

## Testing and development

Run compilation and all automated tests with the interpreter used to install
Annostat:

```bash
python -m compileall -q annostat tests
python -m unittest discover -v
```

For a final end-to-end check, analyze one supplied genome and confirm that
`results/report.html`, `results/summary.json`, the tables, FASTA exports, and SVG
plots are created:

```bash
annostat analyze \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

### Example terminal output

```bash
annostat analyze \
  -f data/GCF_000007145.1.fna \
  -g data/GCF_000007145.1.gff3 \
  -o results \
  --table-format tsv
```

Example output:

```text
Annostat 1.1.0 | bacterial genome annotation analysis
FASTA: /path/to/annostat/data/GCF_000007145.1.fna
GFF3:  /path/to/annostat/data/GCF_000007145.1.gff3

[1/5] Reading GFF3 annotations and FASTA sequences
[2/5] Streaming, translating, and writing CDS sequences
[3/5] Calculating statistics and writing analysis tables
[4/5] Rendering scientific visualizations
[5/5] Building the offline HTML report

Analysis summary
--------------------------------------------------------
  Genome size                              5,076,188 bp
  Sequences                              1 (1 circular)
  Features                                        4,498
  CDS                                             4,295
  RNA features                                      181
  Genome GC                                      65.07%
  Coding density                                 85.33%
  Hypothetical CDS                         619 (14.41%)
  COG-annotated CDS                      2,881 (67.08%)
  Recognized starts                      4,294 (99.98%)
  QC review                1 warning, 872 informational findings
--------------------------------------------------------
  Output                   /path/to/annostat/results
  Files written            14
  Completed in             0.56 seconds
```

Runtime depends on the input size and system. The supplied bacterial genome
typically completes in under a second on a modern desktop. Use `--profile` when
you need measured stage timings and peak-memory information; profiling itself
adds overhead and should not be used for speed benchmarks.

## Scientific basis

The QC checks follow the internal-consistency emphasis of the
[NCBI Prokaryotic Genome Annotation Standards](https://www.ncbi.nlm.nih.gov/refseq/annotation_prok/standards/)
and the attribute conventions in the
[NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).
FASTA validation accepts the NCBI/IUPAC nucleotide alphabet, and translation
tables and initiator codons follow the
[NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
The report treats overlaps as review candidates because overlapping genes are
common in compact bacterial genomes. Translation uses the bacterial code and
alternative-start behavior described by
[Prodigal](https://pmc.ncbi.nlm.nih.gov/articles/PMC2848648/). The project remains
an annotation-analysis tool rather than a replacement for comprehensive
annotation pipelines such as
[Bakta](https://pmc.ncbi.nlm.nih.gov/articles/PMC8743544/) or
[Prokka](https://pubmed.ncbi.nlm.nih.gov/24642063/).

Comparative profiles are normalized to avoid treating larger genomes as
automatically enriched. Full ortholog and pangenome inference remains the domain
of tools such as [Panaroo](https://pmc.ncbi.nlm.nih.gov/articles/PMC7376924/),
while sequence-level average nucleotide identity should be calculated by a
specialized implementation such as
[FastANI](https://pmc.ncbi.nlm.nih.gov/articles/PMC6269478/).
