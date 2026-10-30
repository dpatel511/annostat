# Annostat

[![PyPI](https://img.shields.io/pypi/v/annostat)](https://pypi.org/project/annostat/)
[![CI](https://github.com/dpatel511/annostat/actions/workflows/ci.yml/badge.svg)](https://github.com/dpatel511/annostat/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dpatel511/annostat/blob/main/LICENSE)

Validate, analyze and compare bacterial annotations from matching **FASTA and
GFF3** files. Outputs include offline HTML reports, JSON/TSV results, CDS FASTA
exports and editable SVG plots. Annostat analyzes existing annotations; it does
not annotate genomes or establish biological correctness.

This README describes **1.1.1**. Check `annostat --version`; package indexes
may provide an earlier release. Install this checkout for unreleased fixes.

## Installation

Python **3.10+** is required. Matplotlib and Biopython install automatically.

```bash
python -m pip install annostat
annostat --version
```

Use an isolated environment: `python3 -m venv .venv`, then
`source .venv/bin/activate` on Linux/WSL/macOS or
`.\.venv\Scripts\Activate.ps1` in Windows PowerShell.

Bioconda is available on Linux, WSL and macOS:

```bash
conda create -n annostat -c conda-forge -c bioconda --strict-channel-priority annostat
conda activate annostat
```

For development or an unreleased fix:

```bash
git clone https://github.com/dpatel511/annostat.git
cd annostat
python -m pip install -e .
```

## Quick start

From a source checkout, use the included example files:

```bash
annostat validate -f data/GCF_000007145.1.fna -g data/GCF_000007145.1.gff3
annostat analyze -f data/GCF_000007145.1.fna -g data/GCF_000007145.1.gff3 -o results
```

Open `results/report.html`. Use `summary.json` for machine-readable metrics,
`tables/` for counts and QC findings, `sequences/` for CDS nucleotide/protein
FASTA files, and `plots/` for SVG figures. `--table-format tsv` selects a TSV
feature overview; the default is CSV. `--profile` adds peak Python-memory
measurement and stage timings.

## Choose the comparison mode

| Mode | Input | Result |
|---|---|---|
| `compare` | Two or more annotations | Aggregate counts, normalized profiles and pairwise metric differences |
| `compare --gene-level` | Exactly two annotations of identical genome sequences and IDs | Also matches individual CDS models and compares predicted proteins |

Aggregate comparison:

```bash
annostat compare \
  --genome A genome-a.fna annotation-a.gff3 \
  --genome B genome-b.fna annotation-b.gff3 -o comparison
```

CDS-level comparison on the same unchanged assembly:

```bash
annostat compare --gene-level \
  --genome RefSeq genome.fna refseq.gff3 \
  --genome Bakta genome.fna bakta.gff3 -o gene-comparison
```

Open `gene-comparison/comparison.html` and inspect a CDS to connect its coordinates,
translation context and alignment to the classification. `comparison.txt` is the
text summary; all relationships remain in `tables/gene_comparisons.tsv` and JSON.

**Matching uses CDS structure, strand and frame, not gene names.** Exact unique
matches anchor the comparison, followed by mutually unique start/stop alternatives.
Other overlaps remain unresolved. Equal first-codon letters can occur at different
positions: TTG can initiate a shorter model but encode L internally in a longer
model. Reports distinguish each model's first-codon positions and first predicted
amino acid from aligned shared-codon differences.

See the [gene-comparison reference](https://github.com/dpatel511/annostat/blob/main/docs/gene-comparison.md) for rules, columns,
formulas and limitations, or the [demonstration](https://github.com/dpatel511/annostat/blob/main/docs/gene-level-demo.md) for a
short presentation using real annotation files.

## Other workflows

Filters create additional exports and do not change the full analysis. All enabled
criteria must be satisfied:

```bash
annostat analyze -f genome.fna -g annotation.gff3 -o results \
  --min-cds-length 300 --max-cds-length 900 --require-cog --exclude-hypothetical
```

Aggregate completed single-annotation analyses; sample directory names must be
unique. This does not aggregate gene-comparison tables:

```bash
annostat summarize results/sample-a results/sample-b -o cohort
```

Downloads require the optional official
[NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
and internet access. Local analysis works offline and does not transmit inputs:

```bash
annostat fetch GCF_000007145.1 -o ncbi-data
annostat compare --genome Sample sample.fna sample.gff3 \
  --reference GCF_000007145.1 -o reference-comparison
```

Run `annostat COMMAND --help` for options. `inspect` is an alias for `analyze`;
the older option-only invocation is deprecated. `python -m annostat` is equivalent
to the installed command.

## Scientific scope

- FASTA IDs are the first whitespace-delimited header tokens and must match GFF3
  IDs. Coordinates are **1-based inclusive**; strand and phase matter.
- Repeated CDS IDs count once and join multipart rows; the feature overview
  retains individual GFF3 rows. Phase continuity is validated. Circular
  coordinates extending past the sequence length require a `region` with
  `Is_circular=true`. Gene comparison does not align different or rotated assemblies.
  Reconcile coordinates rather than renaming IDs to bypass checks.
- NCBI tables **4, 11 and 25** are supported. A consistent GFF3 `transl_table` is
  used when present; otherwise the default is 11. Use `analyze --genetic-code TABLE`
  or `compare --genetic-code LABEL TABLE` for inputs without that metadata.
  Conflicting overrides are rejected.
- Initiator M is a translation rule, not experimental evidence. Nonzero initial
  phase and 5′ partial markers prevent initiator recoding. Supported contiguous
  codon `transl_except` annotations include standard amino acids, Sec, Pyl, TERM
  and OTHER (X). In the development checkout, unapplied exception entries produce
  an `unapplied_translation_exception` QC warning; affected accepted pairs have
  limited protein evidence. Their exported translations/lengths are provisional.
  Complex exception positions still need external verification.
- GC is **G+C divided by A+C+G+T**; ambiguity-only genomes have unavailable GC.
  Assembly length includes ambiguous bases. Coding density uses the union of CDS
  intervals, counting overlaps once.
- COG coverage uses CDS with source-provided categories divided by all CDS.
  Category profiles use category assignments; a CDS can have several categories.
  Missing categories are unavailable, not biological zero. Gene-symbol and COG-ID
  Jaccard scores compare labels, not orthologs or ANI.
- First-codon counts include partial CDSs and do not prove initiation. Pooled codon
  usage includes complete codons, including terminal stops; incomplete trailing
  bases are excluded. CDS histograms disclose a lower empirical 99th-percentile view cutoff;
  mean and median use all CDS lengths.
- QC findings identify records for review. Unmatched calls do not establish gene
  absence; annotation differences do not establish DNA mutations, altered function,
  statistical significance or which annotation is correct.

## Outputs and reproducibility

| Workflow | Files to start with |
|---|---|
| `analyze` | `report.html`, `summary.json`, `tables/`, `sequences/`, `plots/` |
| `validate` | `validation.json` and `validation.tsv` when findings exist |
| `compare` | `comparison.html`, `comparison.json`, dataset/pairwise TSV tables and SVG plots |
| `compare --gene-level` | Also `comparison.txt`, `tables/gene_comparisons.tsv`, `plots/gene_consequences.svg` |
| `summarize` | `cohort.html`, `cohort.tsv`, `cohort.json` |

Standalone validation writes nothing for clean inputs. Full analysis retains
`validation/validation.json`; its TSV is written only for findings. COG plots are
omitted when category data is unavailable. Reports work offline; the interactive
gene table requires JavaScript, while TSV/JSON/text exports do not. An analysis
summary's `output_files` list records generated files.

JSON records schemas, input SHA-256 hashes and calculation policy. Scientific
fingerprints exclude paths and timings. Missing measurements remain `null`/`NA`,
separate from zero. Cohorts flag mixed versions or calculation policies;
regenerate inputs consistently before comparing their metrics.

## Terminal output and troubleshooting

Progress and diagnostics go to **stderr**, final summaries to **stdout**:

```bash
annostat analyze -f genome.fna -g annotation.gff3 -o results >summary.txt 2>run.log
```

`--quiet` suppresses routine output; errors remain visible. Exit status **0**
means success, **1** means a validation threshold was reached, and **2** means
invalid options or a runtime failure. Validation defaults to `--fail-on error`;
use `--fail-on warning` for stricter CI or `--fail-on never` to retain findings
without failing a pipeline. Its PASS/FAIL label describes structural validity,
independently of that exit policy.

If the command is missing, activate its environment or use `python -m annostat`.
`--gene-level` belongs to `compare` and needs 1.1.0 or newer. For input errors,
inspect `validate` output before changing files. For bug reports, include the
version, command, error and a small non-sensitive example in the
[issue tracker](https://github.com/dpatel511/annostat/issues).

## Development and references

```bash
python -m compileall -q annostat tests
python -m unittest discover -v
```

Translation and coordinates follow the [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
and [NCBI GFF3 documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).
See the [1.1.1 validation record](https://github.com/dpatel511/annostat/blob/main/docs/validation-1.1.1.md) for the tested panel and limits, and the [changelog](https://github.com/dpatel511/annostat/blob/main/CHANGELOG.md) for release changes. Computational consistency
and source checks do not replace independent biological validation.
