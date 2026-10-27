#!/usr/bin/env bash
# Run from the checkout root, using an environment with Annostat dependencies.
set -euo pipefail

if [ "$#" -ne 5 ]; then
    printf 'Usage: bash examples/gene-level-demo.sh BAKTA_FASTA BAKTA_GFF REFSEQ_FASTA REFSEQ_GFF OUTPUT_DIR\n' >&2
    exit 2
fi
for input in "$1" "$2" "$3" "$4"; do
    if [ ! -s "$input" ]; then
        printf 'Missing or empty input: %s\n' "$input" >&2
        exit 2
    fi
done

inputs=(--genome Bakta "$1" "$2" --genome RefSeq "$3" "$4")
python -m annostat --version
printf '\nAggregate comparison (gene-level toggle off)\n'
python -m annostat compare "${inputs[@]}" -o "$5/aggregate"
printf '\nCDS comparison (gene-level toggle on)\n'
python -m annostat compare "${inputs[@]}" --gene-level -o "$5/gene-level"
python - "$5" <<'PY'
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
aggregate = json.loads((root / 'aggregate/comparison.json').read_text())
gene = json.loads((root / 'gene-level/comparison.json').read_text())['gene_level_comparison']
assert 'gene_level_comparison' not in aggregate
assert sum(gene['consequence_counts'].values()) == gene['status_counts'].get('matched', 0)
print('\nResolved CDS pairs:', gene['status_counts'].get('matched', 0))
print('Protein evidence:', json.dumps(gene['consequence_counts'], sort_keys=True))
print('All relationship rows:', len(gene['relationships']))
print('Genomic coordinates: 1-based inclusive; amino-acid deltas: B minus A.')
print('Experimental support is not assessed by this comparison.')
for path in ('aggregate/comparison.html', 'gene-level/comparison.html', 'gene-level/comparison.txt'):
    print(root.joinpath(path).resolve())
PY
