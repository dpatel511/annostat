"""Check gene-level CLI/report contracts using an installed package, without network access."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
from urllib.parse import quote
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    checkout = Path(__file__).resolve().parents[1]
    location = Path(subprocess.check_output([sys.executable, "-c", "import annostat; print(annostat.__file__)"], text=True).strip()).resolve()
    assert checkout not in location.parents, f"Expected installed package, imported {location}"
    fasta = root / "genome.fna"
    fasta.write_text(">chr\nATGTTGAAATAANNNATGAAAGGGTAA\n")
    for label, start in (("A", 1), ("B", 4)):
        annotation = ";transl_except=" + quote("(pos:join(19..20,21),aa:Leu)", safe="") if label == "A" else ""
        (root / f"{label}.gff3").write_text(
            "##gff-version 3\n" +
            f"chr\ttest\tCDS\t{start}\t12\t.\t+\t0\tID=start_{label}\n" +
            f"chr\ttest\tCDS\t16\t27\t.\t+\t0\tID=exception_{label}{annotation}\n")
    output = root / "comparison"
    command = [sys.executable, "-m", "annostat", "compare"]
    for label in ("A", "B"):
        command += ["--genome", label, str(fasta), str(root / f"{label}.gff3")]
    command += ["-o", str(output)]
    run = subprocess.run(command + ["--gene-level"], capture_output=True, text=True, check=True)
    assert "Accepted CDS pairs" in run.stdout and "gene-level" in run.stderr
    data = json.loads((output / "comparison.json").read_text())
    rows = data["gene_level_comparison"]["relationships"]
    assert len(rows) == 2
    start = next(r for r in rows if r["id_a"] == "start_A")
    assert start["protein_consequence"] == "initiator_recoding"
    assert "A residue 2=L; B residue 1=M" in start["consequence_note"]
    assert start["first_codon_positions_a"] == "1,2,3"
    assert start["first_codon_positions_b"] == "4,5,6"
    limited = next(r for r in rows if r["id_a"] == "exception_A")
    assert limited["protein_identical"] is True
    assert limited["protein_consequence"] == "limited" and limited["alignment_identity_percent"] is None
    assert limited["unapplied_transl_except_a"] == "(pos:join(19..20,21),aa:Leu)"
    with (output / "tables/gene_comparisons.tsv").open(newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["id_a"] for r in table} == {r["id_a"] for r in rows}
    assert next(r for r in table if r["id_a"] == "exception_A")["unapplied_transl_except_a"] == limited["unapplied_transl_except_a"]
    html = (output / "comparison.html").read_text()
    text = (output / "comparison.txt").read_text()
    assert "Unapplied source transl_except" in html and "Unapplied source transl_except" in text
    assert "unapplied_transl_except" in html and "A unapplied transl_except:" in text
    ET.parse(output / "plots/gene_consequences.svg")
    analysis = root / "analysis"
    subprocess.run([sys.executable, "-m", "annostat", "analyze", "-f", str(fasta), "-g", str(root / "A.gff3"), "-o", str(analysis), "--quiet"], capture_output=True, text=True, check=True)
    with (analysis / "tables/annotation_issues.csv").open(newline="") as handle:
        issues = list(csv.DictReader(handle))
    issue = next(r for r in issues if r["issue_type"] == "unapplied_translation_exception")
    assert issue["feature_id"] == "exception_A" and issue["severity"] == "warning"
    summary = json.loads((analysis / "summary.json").read_text())
    assert summary["quality_control"]["issue_counts"]["unapplied_translation_exception"] == 1
    assert "Unapplied Translation Exception" in (analysis / "report.html").read_text()
    quiet = subprocess.run(command + ["--gene-level", "--quiet"], capture_output=True, text=True, check=True)
    assert not quiet.stdout and not quiet.stderr
    aggregate = subprocess.run(command, capture_output=True, text=True, check=True)
    assert "gene_level_comparison" not in json.loads((output / "comparison.json").read_text())
    assert not (output / "comparison.txt").exists()
    assert not (output / "tables/gene_comparisons.tsv").exists()
    subprocess.run(command + ["--gene-level", "--quiet"], capture_output=True, text=True, check=True)
    print(f"PASS installed gene comparison, exception QC, exports, quiet and aggregate cleanup: {location}")


if __name__ == "__main__":
    main()
