# Validation of the 1.1.1 changes

Audit date: 2026-10-28. This records computational checks of the local
`fix/scientific-report-audit-1.1.1` branch, not a published release or a claim
that either annotation is biologically correct.

## Public-data panel

Each row compares frozen Bakta and RefSeq annotations on matching sequence IDs
and identical DNA. Names are taken from the source RefSeq GenBank records.

| Assembly | Source organism | Accepted CDS pairs |
|---|---|---:|
| GCF_000009045.1 | Bacillus subtilis subsp. subtilis str. 168 | 4,138 |
| GCF_000829055.1 | Lacticaseibacillus casei DSM 20011 = JCM 1134 = ATCC 393 | 2,697 |
| GCF_000012745.1 | Thiobacillus denitrificans ATCC 25259 | 2,798 |
| GCF_000218265.1 | Sinorhizobium meliloti SM11 | 6,579 |
| GCF_012044595.1 | Streptomyces antibioticus | 7,367 |
| GCF_000067165.1 | Sorangium cellulosum So ce56 | 9,370 |

The panel contains 35,017 relationship rows and 32,949 accepted pairs. Rows are
not distinct genes. Additional comparisons use the historical HOMD Prokka
annotation of the Bacillus assembly, including a three-way aggregate report.
All 4,214 Prokka published proteins were checked against its FAA and independent
GenBank-based translation.

Sources: [Bakta benchmark](https://zenodo.org/records/5253552),
[NCBI RefSeq assemblies](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/), and
[HOMD Prokka V11.02](https://www.homd.org/ftp/genomes/PROKKA/V11.02/).
Source versions are frozen locally; downloading a newer annotation can change
counts. SHA-256 hashes and byte sizes were recorded for 52 source files.

## Checks and results

- **139 unit tests passed on Python 3.10.22 and 3.13.** Regression tests cover all
  64 canonical codons under tables 4, 11 and 25; both strands; initial phases;
  standard amino-acid exceptions; terminal versus internal stops; and circular
  first-codon positions.
- **102 CLI runs passed:** validation, analysis, filters, quiet mode, forward and
  reverse comparisons, aggregate and gene-level reports, and a 12-sample cohort.
  Small fixtures cover missing versus zero GC, no CDSs, translation-table
  reassignment, unknown residues, internal stops, multipart CDSs and circular
  origins. Expected validation failures were checked for the correct exit status.
- Independent checks use source GenBank coordinates/proteins, FAA, gffutils and
  NumPy without importing Annostat calculations. They check CDS correspondence,
  source sequences, interval-union coding density, GC denominators, filters,
  codon counts, protein alignments and consequence counts.
- All generated reports were checked for finite JSON numbers, TSV/JSON/HTML/text
  agreement, valid local links and offline assets. Profile distances were
  recomputed with 60-digit Decimal arithmetic; pairwise deltas, percentages,
  Jaccard values and cohort source joins were cross-checked.
- **89 HTML reports passed at widths 1440 and 430 pixels.** Interactive row
  selection/filtering, the equal-TTG explanation and safe handling of hostile
  annotation text passed without JavaScript errors or external asset requests.
  **89 distinct SVGs passed browser text-bound checks**; the audit generated
  187 SVG files including repeated runs. Plot-artist tests check plotted values,
  denominators, histogram cutoffs, missing values and heatmap contrast.
- Wheel/source builds and metadata checks passed. A fresh installed wheel on
  Python 3.10 imported outside the checkout and reproduced the complete
  Bacillus gene-comparison JSON from the source run.

Coordinates in reports are **1-based inclusive**, with first-codon bases listed
in translation order. Signed comparison deltas are **B minus A**. Reverse-input
checks preserve correspondence and reverse directional measurements.

## The equal-TTG example

Bakta `KBPOGL_00455` starts at 80844; RefSeq `cds-NP_387953.1` starts at
80802. Each model's first triplet is TTG and its first predicted residue is M.
At the shared triplet 80844,80845,80846, Bakta residue 1 is M under the initiator
rule, whereas RefSeq residue 15 is L under the elongation assignment. Equal
first-triplet strings therefore do not imply equal positions in the proteins.

The report now exposes both codon locations, both first residues and translation
contexts. This case is an interpretation issue; separate tests exposed actual
incomplete-start and translation-exception bugs, which were corrected.

## Scope and limits

Published FAA agreement verifies consistency with source predictions. Pseudogene
CDSs without published proteins are checked as computed translations, separately
from published-protein evidence. Neither check supplies experimental validation.

This panel covers bacterial annotations and the exercised genetic codes and
exception formats. Complex unsupported `transl_except` locations, other organism
groups, every possible malformed annotation and journal-specific figure
requirements are not certified by these tests. Alternative starts need suitable
independent initiation or N-terminal evidence before biological conclusions.
Unmatched calls do not establish gene absence.

## Repeating checks

From an environment with the checkout installed:

```bash
python -m unittest discover -v
python -m compileall -q annostat tests
bash -n examples/gene-level-demo.sh
```

See the [demo](gene-level-demo.md) for the frozen Bacillus inputs and commands.
This audit's full command manifests, source checksums, independent checkers,
logs, reports and screenshots are retained locally under
`tmp/release-1.1.1-audit/`; these large audit resources are ignored by Git and
are not included in the package. They include `source-checksums.json`,
`source-audit-summary.json`, `all-output-verification.json` and
`browser/verification.json`.
