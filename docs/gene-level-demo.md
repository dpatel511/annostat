# Demonstrating gene-level annotation comparison

Use this to show how `compare --gene-level` helps review two annotations of the
same unchanged genome. The toggle adds individual CDS correspondence, translated
sequence differences, source context, and an interactive evidence view.

All genomic coordinates below are **1-based inclusive**. Amino-acid length and
terminal deltas are **B minus A**. This demo measures annotation differences; it
does not select a biologically correct annotation or predict functional impact.

## Prepare before the meeting

Run the current checkout so an older installed release cannot hide the new report.
From WSL, activate an environment containing this checkout's dependencies:

```bash
cd /mnt/e/Projects/annostat
source tmp/external-benchmark/wsl-venv/bin/activate
python -m annostat --version
git branch --show-current
```

The prepared environment and downloaded datasets above are local audit resources,
not files distributed with the repository. For a fresh checkout, create and
activate a virtual environment and run `python -m pip install -e .`, then supply
your own paired inputs to the demo script.

This workspace already has two independently generated annotations of
**Bacillus subtilis subsp. subtilis 168, GCF_000009045.1**: the published Bakta
benchmark annotation and the downloaded RefSeq annotation. Their FASTA sequence IDs
and bases have been cross-checked. Byte hashes of the two FASTA files differ;
sequence equality is checked after parsing, allowing different descriptions and
line wrapping. Neither annotation is treated as experimental truth.

```bash
BAKTA=tmp/external-benchmark/raw/GCF_000009045.1-Bacillus-subtilis-subsp.-subtilis/GCF_000009045.1
REFSEQ=tmp/external-benchmark/refseq-oracle/GCF_000009045.1_ASM904v1
DEMO=tmp/professor-demo

bash examples/gene-level-demo.sh \
  "$BAKTA.fna" "$BAKTA.gff3" \
  "${REFSEQ}_genomic.fna" "${REFSEQ}_genomic.gff" \
  "$DEMO"
```

Prepare the output before presenting. The report works offline, so no download or
annotation run is needed during the meeting. Re-running the commands updates only
the named output folders; keep outputs you want to retain in a separate folder.

## Show the toggle explicitly

These are the two commands used by the script. Use separate output directories so
both reports remain available side by side:

```bash
python -m annostat compare \
  --genome Bakta "$BAKTA.fna" "$BAKTA.gff3" \
  --genome RefSeq "${REFSEQ}_genomic.fna" "${REFSEQ}_genomic.gff" \
  -o "$DEMO/aggregate"

python -m annostat compare \
  --genome Bakta "$BAKTA.fna" "$BAKTA.gff3" \
  --genome RefSeq "${REFSEQ}_genomic.fna" "${REFSEQ}_genomic.gff" \
  --gene-level -o "$DEMO/gene-level"

explorer.exe "$(wslpath -w "$DEMO/aggregate/comparison.html")"
explorer.exe "$(wslpath -w "$DEMO/gene-level/comparison.html")"
```

The aggregate report summarizes annotation profiles. In the gene-level report,
start with **Overview**: the main finding, accepted-pair counts, and per-source
pairing coverage. Then use **Differences** to separate structural boundaries from
protein evidence. **View rows** opens a category in the evidence table and clears
older searches/filters. Use **Inspect** to explain individual differences.

The structural agreement formula uses all CDS models across both sources. Protein
category percentages use accepted pairs, so these percentages have a different
denominator. The protein-length summary reports direction, median absolute change,
and signed range only for interpretable protein differences; it is descriptive,
not a functional ranking or a significance test.

## Four examples to show

Clear the search field and both filters between examples. Search by the exact
CDS ID below if a gene symbol gives more than one hit.

| Search | Finding | What to explain |
|---|---|---|
| `KBPOGL_00010` (`dnaA`) | Exact structure, identical protein; 446 aa in both | An unchanged reference example: coordinates 410–1750, no measured difference. |
| `KBPOGL_00040` (`gyrB`) | Start shift; A 4861–6783, B 4867–6783; 640 versus 638 aa | A contains two extra N-terminal residues, `MA`; all 638 shared residues are identical. The report localizes the disagreement without deciding which start is used in vivo. |
| `KBPOGL_00300` (`prtG`) | Start shift; A 51662–52552, B 51680–52552; 296 versus 290 aa | A has six extra residues, `MCGVRN`; 289 shared residues are identical and one mismatch follows initiation-methionine recoding. A starts at ATG; B starts at GTG. This is not a DNA mutation. |
| `cds-NP_387980.1` (`rpmGB`) | RefSeq-only unmatched CDS at 117349–117498; 49 aa | There is no accepted overlapping Bakta counterpart. That does not establish biological gene absence. Protein-consequence measurements remain unavailable. |

For `gyrB` and `prtG`, show the interpretation, boundary metrics, locus SVG,
changed-sequence window, and source CDS models. Expand **Full protein alignment**
and **Codon-mapped alignment measurements** to inspect all counts; gray highlights
mark codons present in only one model, and amber highlights mark non-gap mismatches.
Try **Download locus SVG**. For `rpmGB`, show how unavailable evidence differs from
a measured zero. Additional source attributes remain under **Source annotations
and database references**.

These are the expected totals for the frozen local input files:

| Measure | Value |
|---|---:|
| Accepted one-to-one CDS pairs | 4,138 |
| Identical translated protein | 3,753 |
| Terminal-only differences | 253 |
| Initiation-context differences | 132 |
| Other translated-sequence differences | 0 |
| Limited protein evidence among accepted pairs | 0 |
| Total relationship rows, including unresolved correspondence | 4,320 |
| Unmatched calls: Bakta / RefSeq | 82 / 95 |
| Candidate / ambiguous relationship rows | 3 / 2 |

Relationship rows are not distinct genes. Consequence counts use accepted pairs
only. Changed inputs or annotation releases can produce different totals.

## Show the text and machine-readable evidence

```bash
cat "$DEMO/gene-level/comparison.txt"

python - "$DEMO/gene-level/tables/gene_comparisons.tsv" <<'PY'
import csv
import sys

targets = {'KBPOGL_00010', 'KBPOGL_00040', 'KBPOGL_00300', 'cds-NP_387980.1'}
fields = ('gene_a', 'gene_b', 'id_a', 'id_b', 'relationship', 'protein_consequence',
          'protein_length_delta_aa', 'initiator_mismatch_count')
print('\t'.join(fields))
with open(sys.argv[1], newline='', encoding='utf-8') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        if targets.intersection((row['id_a'], row['id_b'])):
            print('\t'.join(row[field] or 'NA' for field in fields))
PY
```

The text shortlist includes examples from limited evidence and unresolved matches,
not just the largest measured length differences. It shows at most three examples
per group; the HTML, TSV, and JSON retain all relationships.

## Optional: show limited evidence on a second assembly

The Bacillus accepted pairs have no limited translations. The local six-assembly
audit also includes a Lactobacillus casei pair with 54 limited accepted pairs:

```bash
CASE=tmp/external-benchmark/panel-raw/GCF_000829055.1
python -m annostat compare \
  --genome Bakta "$CASE/bakta/GCF_000829055.1.fna" "$CASE/bakta/GCF_000829055.1.gff3" \
  --genome RefSeq "$CASE/GCF_000829055.1_ASM82905v1_genomic.fna" "$CASE/GCF_000829055.1_ASM82905v1_genomic.gff" \
  --gene-level -o "$DEMO/limited-evidence"
explorer.exe "$(wslpath -w "$DEMO/limited-evidence/comparison.html")"
```

Choose **Limited protein evidence** and search `HPIBKM_00040`. Its counterpart is
`cds-LBCZ_RS00035`: the predictions have 104 and 126 aa, but unknown residues or
internal stops prevent reliable residue-level consequence assignment. RefSeq
provides `pseudo=true`; this is source context, not experimental proof. The length
difference is measurable, while identity and residue-mismatch counts remain
unavailable. This shows why uncertainty should remain visible in a useful report.

## Optional: reverse the comparison

```bash
python -m annostat compare \
  --genome RefSeq "${REFSEQ}_genomic.fna" "${REFSEQ}_genomic.gff" \
  --genome Bakta "$BAKTA.fna" "$BAKTA.gff3" \
  --gene-level -o "$DEMO/reversed"
```

Search `KBPOGL_00040` in the reversed report: the protein-length delta changes
from −2 to +2 aa because the source order changes; the underlying disagreement
and consequence counts stay the same.

## Explain the biological follow-up

The report answers **where annotations disagree and how their predicted proteins
differ**. It does not assess experimental support. For start-site differences,
independent initiation-focused ribosome profiling or diagnostic N-terminal peptide
evidence can support one or both alternatives. Peptides shared by both models do
not distinguish their starts; non-detection alone does not establish absence.

For a published example of investigating start-site changes using ribosome
profiling, N-terminal proteomics, and reporter evidence, see
[Nakahigashi et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4909307/).
This paper motivates a separate biological benchmark; its evidence has not been
mapped onto or used to validate the Bacillus examples above. An evidence benchmark
must verify strain, assembly, coordinates, detection coverage, and the source
study's identification/error-control thresholds before drawing conclusions.
