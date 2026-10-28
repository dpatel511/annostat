# Gene-level CDS comparison

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

## Read the evidence

For a runnable WSL demonstration, verified real-data examples, and an explanation
of what to show, see the [gene-level walkthrough](gene-level-demo.md).

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

Gene-level JSON schema **1.3** records these fields, `consequence_counts`,
`unmatched_cds_counts`, and method details. Matching/status rules and existing columns are retained; first-codon context
fields are documented below. In TSV, null values are blank; in JSON they are `null`.
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


## First codons and amino acids

First codons belong to each model's own coding boundary. Equal DNA triplets are
not necessarily the same genomic occurrence. In the frozen Bacillus example,
RefSeq `yacD` and Bakta `SurA` both begin with TTG, at positions 80802 and 80844
respectively (1-based inclusive). At bases 80844,80845,80846, RefSeq residue 15
is L, while Bakta residue 1 is M under the initiator rule. Both models' own first
amino acid is M. This does not establish a DNA substitution.

| Paired fields (`_a` / `_b`) | Meaning |
|---|---|
| `start_codon` | First complete translated DNA triplet after phase adjustment; not proof of initiation |
| `first_codon_positions` | Three genomic bases in translation order, 1-based; may cross the circular origin |
| `first_amino_acid` | First residue actually present in the predicted protein |
| `start_context` | `initiator`, `incomplete_start`, `translation_exception`, `unrecognized_codon`, or `unavailable` |

Initiation-context notes identify the shared genomic codon and the different
1-based protein positions. Nonzero initial phase or a 5′ partial marker prevents
initiator recoding; an explicit first-codon exception takes precedence. The first
codon labels do not describe every codon in the alignment.

The gene schema is 1.3 in 1.1.1; these context fields also appear in the TSV.
Calculation policy `annotation-metrics-1.2` records the scientific corrections.
Analysis/comparison/cohort schemas remain 1.1.
