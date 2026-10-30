# Repeatable public-data comparison

`examples/public-comparison-manifest.json` freezes four uncompressed inputs from
the 2026-10-28 Bacillus subtilis 168 audit (GCF_000009045.1): Bakta and RefSeq
FASTA/GFF3 pairs. Its source record URLs, original member names, byte sizes and
SHA-256 values identify the exact annotation snapshots. The four files total
about 18 MB (17.5 MiB). Data are not distributed in this repository, and the checker never
downloads or edits inputs. Download/extract the source records named in the
manifest, then name the selected files `bakta.fna`, `bakta.gff3`, `refseq.fna` and
`refseq.gff3` in one input directory. RefSeq downloads may require gzip decompression.

From an environment with the intended checkout/wheel installed:

```bash
python examples/check_public_comparison.py \
  --inputs /path/to/frozen-inputs \
  --output tmp/public-comparison-check
```

A changed annotation or wrong file fails its checksum check before analysis;
do not replace hashes or expected results merely to make it pass. Source archives
can change or become unavailable. A new snapshot requires an independently reviewed
manifest and expectations, not an assumption that the old audit still applies.

The checker verifies 4,138 accepted pairs, 4,320 relationship rows, frozen protein
category counts, alignment identity arithmetic, TSV ID joins and the shared TTG
example at 80844,80845,80846 (1-based inclusive). It writes
`public-data-verification.json` with input hashes, package/policy/schema versions
and results. Outputs overwrite the chosen comparison directory as in normal CLI use.
This is a compact frozen regression, not a replacement for the six-assembly audit,
independent source-protein checks or experimental biological validation. It does
not assess all report layouts or other organisms.

CI uses `examples/check_installed_gene_comparison.py` separately: a small offline
synthetic comparison from an isolated installed wheel checks initiation context,
unapplied exceptions, QC, HTML/text/TSV/JSON exports, quiet mode and removal of
stale gene-level outputs on an aggregate rerun. It downloads no public data.

```bash
# Run from outside the repository with the wheel installed in this Python environment.
python /path/to/annostat/examples/check_installed_gene_comparison.py \
  --output /path/to/smoke-results
```
