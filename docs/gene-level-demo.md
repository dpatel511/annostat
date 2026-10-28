# Demonstrate gene-level comparison

Show two annotations of the **same unchanged genome**. Coordinates below are
**1-based inclusive**; protein-length deltas are **B minus A**. The comparison
measures annotation differences, not experimental correctness or functional impact.

## Prepare the reports

Activate an environment with this checkout installed (`python -m pip install -e .`).
The demo script takes four input files and an output directory:

```bash
python -m annostat --version
bash examples/gene-level-demo.sh bakta.fna bakta.gff3 refseq.fna refseq.gff3 demo
```

It runs `compare` once without and once with `--gene-level`, retaining both
reports under `demo/aggregate/` and `demo/gene-level/`. Prepare them before the
meeting; the HTML works offline. Reusing an output directory replaces its reports.

The following frozen local inputs were independently annotated as **Bacillus
subtilis 168, GCF_000009045.1**. They are audit resources, not distributed files:

```bash
BAKTA=tmp/external-benchmark/raw/GCF_000009045.1-Bacillus-subtilis-subsp.-subtilis/GCF_000009045.1
REFSEQ=tmp/external-benchmark/refseq-oracle/GCF_000009045.1_ASM904v1
bash examples/gene-level-demo.sh \
  "$BAKTA.fna" "$BAKTA.gff3" \
  "${REFSEQ}_genomic.fna" "${REFSEQ}_genomic.gff" demo
```

Sources: [published Bakta benchmark](https://zenodo.org/records/5253552) and
[RefSeq assembly files](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/).
Their parsed sequence IDs and bases agree; FASTA formatting can differ. New source
annotation releases can change the totals below.

## Five-minute walkthrough

1. Open `demo/aggregate/comparison.html`: it compares profiles, without individual
   CDS matching. Then open `demo/gene-level/comparison.html`.
2. Read **Overview** and **Differences**. Protein-category percentages use accepted
   pairs; all-model structural agreement has a different denominator.
3. Use **View rows**, search a CDS ID below, then **Inspect**. Clear old searches
   and filters between examples.
4. Connect the interpretation to the locus, source models, first-codon positions
   and protein alignment. Expand measurements/full alignment when needed.
5. Show `demo/gene-level/comparison.txt` and the complete TSV/JSON exports.

| Search ID | Finding to explain |
|---|---|
| `KBPOGL_00010` (`dnaA`) | Same CDS structure and 446-aa protein. Agreement does not independently prove correctness. |
| `KBPOGL_00040` (`gyrB`) | Bakta 640 aa, RefSeq 638 aa: delta −2, with extra `MA` in Bakta; shared residues match. |
| `KBPOGL_00300` (`prtG`) | Bakta 296 aa, RefSeq 290 aa: delta −6; one shared-codon mismatch follows initiation context. |
| `KBPOGL_00455` (Bakta `SurA`, RefSeq `yacD`) | Both first codons are TTG, at 80844 and 80802. At bases 80844,80845,80846, Bakta residue 1 is M and RefSeq residue 15 is L. These are different initiation contexts on identical DNA. |
| `cds-NP_387980.1` (`rpmGB`) | RefSeq-only call: no accepted Bakta counterpart, rather than proof of biological absence. |

The first-codon table describes **each model's own beginning**. It does not say
that every aligned mismatch occurs at those two first codons. The initiation
note names the shared genomic codon and both protein positions.

Expected totals for these frozen inputs:

| Measure | Count |
|---|---:|
| Accepted pairs | 4,138 |
| Identical proteins | 3,753 |
| Terminal-only differences | 253 |
| Initiation-context differences | 132 |
| Other translation / limited evidence among accepted pairs | 0 / 0 |
| All relationship rows | 4,320 |
| Unmatched Bakta / RefSeq calls | 82 / 95 |
| Candidate / ambiguous rows | 3 / 2 |

Relationship rows are not distinct genes. Reversing input order reverses signed
deltas, while the underlying correspondence and consequence counts remain.

## Limits and biological follow-up

Limited evidence keeps unknown/internal-stop translations and uncertain codon
mapping visible; unavailable identity is not 0%. Source pseudogene/partial flags
are context, not experimental confirmation. See the
[reference](gene-comparison.md) for definitions, fields and formulas.

Independent initiation-focused ribosome profiling or diagnostic N-terminal
peptides can help investigate alternative starts. Shared peptides cannot
separate them, and non-detection alone is insufficient. An experimental example is
[Nakahigashi et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4909307/);
its evidence has **not** been used to validate these Bacillus predictions.
