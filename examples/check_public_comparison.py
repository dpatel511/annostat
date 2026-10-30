"""Offline regression of a frozen public Bakta/RefSeq comparison; see docs/public-data-check.md."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True, help="directory containing the four manifest-named uncompressed files")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(Path(__file__).with_name("public-comparison-manifest.json").read_text())
    inputs = args.inputs.resolve()
    for name, info in manifest["files"].items():
        path = inputs / name
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if path.stat().st_size != info["bytes"] or digest.hexdigest() != info["sha256"]:
            parser.error(f"Frozen source mismatch: {path}; do not update expected results to conceal a changed annotation")
    command = [sys.executable, "-m", "annostat", "compare", "--gene-level"]
    for name, label in (("bakta", "Bakta"), ("refseq", "RefSeq")):
        command += ["--genome", label, str(inputs / f"{name}.fna"), str(inputs / f"{name}.gff3")]
    output = args.output.resolve()
    subprocess.run(command + ["-o", str(output)], check=True)
    data = json.loads((output / "comparison.json").read_text())
    gene = data["gene_level_comparison"]
    rows = gene["relationships"]
    expected = manifest["expected"]
    assert len(rows) == expected["relationship_rows"]
    matched = [r for r in rows if r["status"] == "matched"]
    assert len(matched) == expected["accepted_pairs"]
    assert dict(Counter(r["protein_consequence"] for r in matched)) == expected["consequence_counts"]
    for row in matched:
        if row["alignment_identity_percent"] is not None:
            columns = row["identical_aa"] + row["substitution_count"] + row["gap_aa"]
            assert abs(row["alignment_identity_percent"] - 100 * row["identical_aa"] / columns) < 1e-9
    example = next(r for r in rows if r["id_a"] == "KBPOGL_00455")
    assert example["id_b"] == "cds-NP_387953.1"
    assert example["first_codon_positions_a"] == "80844,80845,80846"
    assert example["first_codon_positions_b"] == "80802,80803,80804"
    assert example["first_amino_acid_a"] == example["first_amino_acid_b"] == "M"
    assert "A residue 1=M; B residue 15=L" in example["consequence_note"]
    with (output / "tables/gene_comparisons.tsv").open(newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert [(r["id_a"],r["id_b"],r["status"]) for r in table] == [(r["id_a"] or "",r["id_b"] or "",r["status"]) for r in rows]
    evidence = {"manifest":manifest,"annostat_version":data["annostat_version"],
                "calculation_policy":data["calculation_policy"],"gene_schema":gene["schema_version"],
                "accepted_pairs":len(matched),"relationship_rows":len(rows),
                "limitation":"Frozen regression/source consistency, not experimental correctness or comprehensive biological validation."}
    (output / "public-data-verification.json").write_text(json.dumps(evidence,indent=2)+"\n")
    print("PASS frozen public-source hashes, pairing totals, consequence counts, identity arithmetic, TTG context and TSV IDs")


if __name__ == "__main__":
    main()
