"""Independent checks of Annostat analyses/comparisons of spliced gene/mRNA/CDS annotations.

Imports no Annostat calculations. GFF3 coordinates are 1-based inclusive and are
converted explicitly to Biopython's 0-based half-open locations. A CDS model is
one CDS ID; its rows are ordered in transcription direction and the first row's
phase trims the joined sequence.

Examples:
  python examples/check_fungal_comparison.py --genetic-code 1 \
    --analysis results/a genome.fa a.gff3 --source-proteins a.proteins.fa \
    --comparison results/self A genome.fa a.gff3 B genome.fa a.gff3 --expect-self
"""
import argparse
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

failures: list[str] = []


def check(condition, message):
    if not condition:
        failures.append(message)
    return condition


def read_fasta(path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def read_gff(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            attributes = {}
            for item in fields[8].strip().strip(";").split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key.strip()] = unquote(value)
            rows.append(dict(seqid=fields[0], type=fields[2], start=int(fields[3]), end=int(fields[4]),
                             strand=fields[6], phase=None if fields[7] == "." else int(fields[7]),
                             attributes=attributes))
    return rows


def models(rows):
    """Group CDS rows by ID into ordered blocks; also return hierarchy counts."""
    groups = defaultdict(list)
    for row in rows:
        if row["type"] == "CDS":
            groups[row["attributes"].get("ID") or f"row{id(row)}"].append(row)
    result = {}
    for cds_id, parts in groups.items():
        minus = parts[0]["strand"] == "-"
        ordered = sorted(parts, key=lambda part: part["start"], reverse=minus)
        result[cds_id] = ordered
    unique = lambda kind: len({row["attributes"].get("ID") for row in rows if row["type"] == kind})
    counts = dict(gene=unique("gene"), mRNA=unique("mRNA"), cds_models=len(result),
                  cds_rows=sum(row["type"] == "CDS" for row in rows))
    return result, counts


def five_prime_partial(first):
    attributes = first["attributes"]
    boundary = "end_range" if first["strand"] == "-" else "start_range"
    generic = attributes.get("partial", "").lower() in {"1", "true", "yes"} and not ({"start_range", "end_range"} & attributes.keys())
    return bool(first["phase"]) or boundary in attributes or generic


def coding_sequence(ordered, genome):
    sign = -1 if ordered[0]["strand"] == "-" else 1
    parts = [SimpleLocation(row["start"] - 1, row["end"], strand=sign) for row in ordered]
    location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
    return str(SeqFeature(location).extract(Seq(genome[ordered[0]["seqid"]])))[ordered[0]["phase"] or 0:]


def cds_table(ordered, default, region_tables=None):
    """A CDS's own transl_table, else its sequence region's, else the selected default."""
    own = {int(row["attributes"]["transl_table"]) for row in ordered if "transl_table" in row["attributes"]}
    return own.pop() if own else (region_tables or {}).get(ordered[0]["seqid"], default)


def region_tables(rows):
    return {row["seqid"]: int(row["attributes"]["transl_table"]) for row in rows
            if row["type"] == "region" and "transl_table" in row["attributes"]}


def translate(ordered, genome, default, regions=None):
    """Whole codons only; codons with non-ACGT bases become X, as Annostat documents."""
    table = cds_table(ordered, default, regions)
    coding = coding_sequence(ordered, genome)
    whole = coding[: len(coding) - len(coding) % 3]
    masked = "".join(codon if set(codon) <= set("ACGT") else "NNN"
                     for codon in (whole[i:i + 3] for i in range(0, len(whole), 3)))
    protein = str(Seq(masked).translate(table=table))
    if not five_prime_partial(ordered[0]) and whole[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
        protein = "M" + protein[1:]
    return protein.removesuffix("*")


def has_exception(ordered):
    return any("transl_except" in row["attributes"] for row in ordered)


def signature(ordered):
    first = ordered[0]
    return (first["seqid"], first["strand"], first["phase"] or 0,
            tuple((row["start"], row["end"]) for row in ordered))


def check_analysis(directory, fasta, gff, table, source_proteins=None):
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    genome = read_fasta(fasta)
    rows = read_gff(gff)
    cds, counts = models(rows)
    regions = region_tables(rows)
    name = directory.name
    check(summary["genetic_code"] == table, f"{name}: genetic code {summary['genetic_code']} != {table}")
    tables = Counter(str(cds_table(ordered, table, regions)) for ordered in cds.values())
    check(summary["genetic_code_counts"] == dict(sorted(tables.items(), key=lambda item: int(item[0]))),
          f"{name}: per-CDS table counts {summary['genetic_code_counts']} != {dict(tables)}")
    check(summary["cds_count"] == counts["cds_models"], f"{name}: CDS models {summary['cds_count']} != {counts['cds_models']}")
    check(summary["cds_row_count"] == counts["cds_rows"], f"{name}: CDS rows {summary['cds_row_count']} != {counts['cds_rows']}")
    for kind in ("gene", "mRNA"):
        check(summary["feature_type_counts"].get(kind, 0) == counts[kind], f"{name}: {kind} count differs")
    exported = read_fasta(directory / "sequences" / "cds_protein.fasta")
    check(set(exported) == set(cds), f"{name}: exported protein IDs differ from CDS IDs")
    # transl_except recoding is checked against the source proteins, not re-implemented here.
    exceptions = {cds_id for cds_id, ordered in cds.items() if has_exception(ordered)}
    mismatches = [cds_id for cds_id, ordered in cds.items() if cds_id in exported and cds_id not in exceptions
                  and exported[cds_id] != translate(ordered, genome, table, regions)]
    check(not mismatches, f"{name}: {len(mismatches)} proteins differ from Biopython, e.g. {mismatches[:3]}")
    source = Counter()
    if source_proteins:
        reference = {key: value.rstrip("*") for key, value in read_fasta(source_proteins).items()}
        for cds_id, ordered in cds.items():
            attributes = ordered[0]["attributes"]
            key = next((k for k in (attributes.get("protein_id"), attributes.get("Parent"), cds_id) if k in reference), None)
            ours, theirs = exported.get(cds_id, ""), reference.get(key)
            trailing = len(coding_sequence(ordered, genome)) % 3
            if theirs is None:
                kind = "no_source_protein"
            elif theirs == ours:
                kind = "identical"
            elif trailing and theirs[:-1] == ours:
                kind = "source_translates_partial_final_codon"  # GenBank convention; Annostat drops trailing bases
            elif len(theirs) == len(ours) and all(a == b or a == "X" for a, b in zip(ours, theirs)):
                kind = "source_resolves_ambiguous_codon"  # Annostat keeps X for non-ACGT codons
            else:
                kind = "different"
                failures.append(f"{name}: {cds_id} differs from source protein {key}")
            source[kind] += 1
    return dict(analysis=name, counts=counts, genetic_code_counts=dict(tables), internal_stop_proteins=sum("*" in p for p in exported.values()),
                biopython_mismatches=len(mismatches), transl_except_models=len(exceptions), source_proteins=dict(source))


def check_comparison(directory, left, right, table, expect_self=False):
    data = json.loads((directory / "comparison.json").read_text(encoding="utf-8"))
    gene = data["gene_level_comparison"]
    name = directory.name
    sides = []
    for label, fasta, gff in (left, right):
        rows = read_gff(gff)
        cds, _ = models(rows)
        sides.append((label, cds, read_fasta(fasta), region_tables(rows)))
    check([d["genetic_code"] for d in data["datasets"]] == [table, table], f"{name}: dataset default genetic codes differ")
    for label, cds, _, _ in sides:
        check(gene["cds_counts"][label] == len(cds), f"{name}: {label} CDS count {gene['cds_counts'][label]} != {len(cds)}")
    signatures = [Counter(signature(ordered) for ordered in cds.values()) for _, cds, _, _ in sides]
    exact = sum(min(count, signatures[1][key]) for key, count in signatures[0].items())
    check(gene["exact_model_count"] == exact, f"{name}: exact models {gene['exact_model_count']} != {exact}")
    proteins = [{cds_id: translate(ordered, genome, table, regions) for cds_id, ordered in cds.items()}
                for _, cds, genome, regions in sides]
    accepted = [row for row in gene["relationships"] if row["status"] == "matched"]
    stopped = [row for row in accepted if "*" in proteins[0][row["id_a"]] or "*" in proteins[1][row["id_b"]]]
    check(all(row["protein_consequence"] == "limited" for row in stopped), f"{name}: internal-stop pair not limited")
    identical = [row for row in accepted if row["protein_consequence"] == "identical"]
    check(all(proteins[0][row["id_a"]] == proteins[1][row["id_b"]] for row in identical), f"{name}: identical pair proteins differ")
    if expect_self:
        duplicated = {key for key, count in signatures[0].items() if count > 1}
        unpaired = {row["id_a"] for row in gene["relationships"] if row["status"] != "matched"}
        by_id = sides[0][1]
        check(all(signature(by_id[cds_id]) in duplicated for cds_id in unpaired),
              f"{name}: self-comparison leaves a non-duplicate structure unpaired")
        check(len(accepted) == len(by_id) - sum(signatures[0][key] for key in duplicated),
              f"{name}: self-comparison accepted pairs do not equal unique structures")
    return dict(comparison=name, cds_counts=gene["cds_counts"], status_counts=gene["status_counts"],
                exact_model_count=exact, accepted=len(accepted), internal_stop_accepted=len(stopped),
                consequence_counts=gene["consequence_counts"])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genetic-code", type=int, default=1,
                        help="table for CDSs without their own or a region transl_table (default 1)")
    parser.add_argument("--analysis", nargs=3, action="append", default=[], metavar=("DIR", "FASTA", "GFF3"))
    parser.add_argument("--source-proteins", action="append", default=[],
                        help="protein FASTA for the matching --analysis (same order); IDs are protein_id, Parent or CDS ID")
    parser.add_argument("--comparison", nargs=7, action="append", default=[],
                        metavar=("DIR", "LABEL_A", "FASTA_A", "GFF3_A", "LABEL_B", "FASTA_B", "GFF3_B"))
    parser.add_argument("--expect-self", action="store_true", help="comparisons are self-comparisons")
    parser.add_argument("--output", type=Path, help="write results JSON here")
    args = parser.parse_args()
    results = []
    for index, (directory, fasta, gff) in enumerate(args.analysis):
        source = args.source_proteins[index] if index < len(args.source_proteins) else None
        results.append(check_analysis(Path(directory), fasta, gff, args.genetic_code, source))
    for directory, la, fa, ga, lb, fb, gb in args.comparison:
        results.append(check_comparison(Path(directory), (la, fa, ga), (lb, fb, gb), args.genetic_code, args.expect_self))
    report = dict(results=results, failures=failures)
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
    print("PASS" if not failures else f"FAIL: {len(failures)} check(s)", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
