"""Run and independently audit a frozen, multi-tool public annotation panel.

Run from the source checkout using the intended Python environment. See
 docs/annotation-panel-testing.md for dependencies, source preparation and limits.
"""
import argparse
import csv
import hashlib
import itertools
import importlib.metadata
import platform
import json
import math
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from pathlib import Path
from statistics import median
from xml.etree import ElementTree as ET

from Bio.Data import CodonTable

from annotation_panel_oracle import check, check_analysis, check_pair, check_figures, fasta, load_models


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def json_file(path):
    return json.loads(path.read_text())


def check_exports(folder, gene):
    with (folder/'tables/gene_comparisons.tsv').open(newline='') as f:
        table = list(csv.DictReader(f, delimiter='\t'))
    check(table == [{k: '' if v is None else str(v) for k, v in r.items()} for r in gene['relationships']], 'all TSV/JSON gene fields differ')
    html = (folder/'comparison.html').read_text()
    payload = json.loads(re.search(r'<script id="gene-data" type="application/json">(.*?)</script>', html, re.S).group(1))
    check(payload == gene, 'HTML scientific payload differs')
    groups = json.loads(re.search(r'<script id="gene-loci" type="application/json">(.*?)</script>', html, re.S).group(1))
    unresolved = {i for i, r in enumerate(gene['relationships']) if r['status'] != 'matched'}
    included = [i for group in groups for i in group['rows']]
    check(set(included) == unresolved and len(included) == len(unresolved), 'group membership omits/duplicates unresolved rows')
    # A model belongs to exactly one connected component; each group is connected.
    owner = {}
    for index, group in enumerate(groups):
        nodes, edges = set(), []
        for i in group['rows']:
            row = gene['relationships'][i]
            members = {(s, row['seqid'], row['strand'], row['id_'+s]) for s in ('a', 'b') if row['id_'+s] is not None}
            nodes |= members
            edges.append(members)
        for node in nodes:
            check(node not in owner, 'CDS is split across review groups')
            owner[node] = index
        visited = set(next(iter(edges)))
        while True:
            expanded = visited | set().union(*(e for e in edges if e & visited))
            if expanded == visited:
                break
            visited = expanded
        check(visited == nodes, 'disconnected review group')
        check(group['cds_counts'] == {s: sum(n[0] == s for n in nodes) for s in ('a', 'b')}, 'distinct CDS/group counts differ')
        statuses = {gene['relationships'][i]['status'] for i in group['rows']}
        kind = 'ambiguous' if 'ambiguous' in statuses else 'candidate' if 'candidate' in statuses else 'unmatched'
        check(group['kind'] == kind, 'review group classification differs')
    text = (folder/'comparison.txt').read_text()
    check(f"Accepted one-to-one CDS pairs: {gene['status_counts'].get('matched', 0):,}" in text, 'text accepted denominator differs')
    check('1-based inclusive' in text and 'B minus A' in text, 'text coordinate/direction disclosure missing')
    for svg in (folder/'plots').glob('*.svg'):
        tree = ET.parse(svg)
        check(tree.find('{http://www.w3.org/2000/svg}title') is not None, 'SVG accessible title missing')
    check_figures(folder, gene=gene, profiles=json_file(folder/'comparison.json')['datasets'])
    return sum(group['kind'] == 'ambiguous' for group in groups)


def check_terminal(path, expected):
    text = path.read_text()
    for label, value in expected.items():
        match = re.search(r'^\s*'+re.escape(label)+r'\s+(.+)$', text, re.M)
        check(match is not None and match.group(1).strip() == value, 'terminal value differs: '+label+' in '+str(path))


def check_analysis_terminal(path, summary):
    s, q = summary, summary['quality_control']
    n = s['cds_count']
    recognized = sum(v for k,v in s['start_codon_counts'].items() if k in CodonTable.unambiguous_dna_by_id[s['genetic_code']].start_codons)
    expected = {'Genome size':f"{s['genome_length']:,} bp", 'CDS':f'{n:,}', 'Features':f"{s['total_features']:,}",
                'RNA features':f"{sum(s['rna_counts'].values()):,}", 'Genome GC':f"{q['genome_gc_percent']:.2f}%",
                'Coding density':f"{q['coding_density_percent']:.2f}%", 'Recognized starts':f'{recognized:,} ({100*recognized/n:.2f}%)',
                'Hypothetical CDS':f"{s['hypothetical_cds_count']:,} ({100*s['hypothetical_cds_count']/n:.2f}%)",
                'COG-annotated CDS':f"{s['cds_with_cog_count']:,} ({100*s['cds_with_cog_count']/n:.2f}%)" if s['cog_data_available'] else 'not available',
                'Files written':str(len(s['output_files']))}
    check_terminal(path,expected)


def check_profile(profile, summary, models):
    for name in ('cds_count', 'genome_length', 'genetic_code', 'calculation_policy'):
        check(profile[name] == summary[name], 'aggregate/analysis profile differs: '+name)
    check(profile['median_cds_length'] == median(len(m['path']) for m in models.values()), 'median CDS length differs')
    for key, value in [('gc_percent', summary['quality_control']['genome_gc_percent']), ('coding_density_percent', summary['quality_control']['coding_density_percent'])]:
        check(math.isclose(profile[key], value), 'aggregate/analysis metric differs: '+key)



def check_pairwise_metrics(document):
    profiles = {p['label']: p for p in document['datasets']}
    for row in document['pairwise_comparisons']:
        a, b = profiles[row['dataset_a']], profiles[row['dataset_b']]
        for key, value in row.items():
            if not key.endswith('_delta') or key.endswith('_percent_delta'):
                continue
            metric = key[:-6]
            expected = None if a[metric] is None or b[metric] is None else b[metric]-a[metric]
            check(value == expected, 'aggregate signed delta differs: '+metric)
            percentage = None if expected is None or not a[metric] else 100*expected/a[metric]
            check(row[metric+'_percent_delta'] == percentage, 'relative delta denominator differs: '+metric)
        for key, metric in [('gene_symbol_jaccard','gene_symbols'), ('cog_id_jaccard','cog_ids')]:
            x, y = set(a[metric]), set(b[metric])
            expected = len(x & y)/len(x | y) if x | y else None
            if metric == 'cog_ids' and not (x and y):
                expected = None
            check(row[key] == expected, 'Jaccard/missing-data handling differs')
        for key, metric in [('start_profile_distance','start_codon_percentages'), ('cog_profile_distance','cog_category_percentages')]:
            x, y = a[metric], b[metric]
            if not (sum(x.values()) and sum(y.values())):
                expected = None
            else:
                labels = sorted(x.keys() | y.keys())
                p = [x.get(k,0)/sum(x.values()) for k in labels]
                q = [y.get(k,0)/sum(y.values()) for k in labels]
                expected = math.sqrt(max(0, sum((v*math.log2(2*v/(v+w)) if v else 0)+(w*math.log2(2*w/(v+w)) if w else 0) for v,w in zip(p,q))/2))
            check(row[key] is None if expected is None else math.isclose(row[key], expected, abs_tol=1e-12), 'profile distance differs')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True, help='prepared panel directory containing input-manifest.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--check-only', action='store_true', help='audit existing outputs without rerunning CLI commands')
    args = parser.parse_args()
    root, output = args.inputs.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json_file(root/'input-manifest.json')
    check(manifest == json_file(Path(__file__).with_name('annotation-panel-manifest.json')), 'prepared panel manifest differs from the frozen repository manifest')
    runs, findings, comparisons, analyses = [], [], [], []
    log_lock = Lock()
    environment = subprocess.run([sys.executable, '-c', 'import annostat,sys; print(sys.version); print(annostat.__version__); print(annostat.__file__)'], capture_output=True, text=True, check=True).stdout
    (output/'environment.txt').write_text(environment)
    repo = Path(__file__).resolve().parents[1]
    revision = subprocess.run(['git','-C',str(repo),'rev-parse','HEAD'],capture_output=True,text=True)
    state = dict(python=sys.version, platform=platform.platform(), source_commit=revision.stdout.strip(),
                 modules={p.relative_to(repo).as_posix(): digest(p) for p in sorted((repo/'annostat').rglob('*'))
                          if p.is_file() and p.suffix in {'.py','.js','.css'}},
                 dependencies={name: importlib.metadata.version(name) for name in ('biopython','matplotlib','numpy','gffutils','pyrodigal')})
    (output/'environment.json').write_text(json.dumps(state,indent=2)+'\n')

    def run(name, arguments, quiet=False, expected=0):
        if args.check_only:
            return
        command = [sys.executable, '-m', 'annostat', *map(str, arguments)]
        before = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True, timeout=600)
        log = output/'logs'/name
        log.parent.mkdir(parents=True, exist_ok=True)
        Path(str(log)+'.stdout').write_text(result.stdout)
        Path(str(log)+'.stderr').write_text(result.stderr)
        record = dict(name=name, command=command, returncode=result.returncode, expected=expected,
                      seconds=time.perf_counter()-before, stdout_bytes=len(result.stdout), stderr_bytes=len(result.stderr))
        with log_lock:
            runs.append(record)
            (output/'cli-runs.json').write_text(json.dumps(runs, indent=2)+'\n')
        check(result.returncode == expected, f'{name}: expected exit {expected}, got {result.returncode}: {result.stderr[-2000:]}')
        check('Traceback' not in result.stderr and 'Traceback' not in result.stdout, name+': traceback')
        if quiet and expected == 0:
            check(not result.stdout and not result.stderr, name+': quiet output is not empty')
        if not quiet and expected == 0:
            check(result.stdout.strip(), name+': missing result summary')
            if arguments[0] in {'analyze', 'compare'}:
                check('Annostat' in result.stderr, name+': missing progress stream')
        print('CLI', name, 'PASS', flush=True)

    for case in manifest['cases']:
        acc, sources = case['accession'], case['tools']
        case_result = dict(accession=acc, organism=case['name'], tools={}, pairs=[])
        models, summaries = {}, {}
        try:
            genome = fasta(root/sources['RefSeq']['fasta'])
            for label, source in sources.items():
                for kind, expected_hash in source['sha256'].items():
                    check(digest(root/source[kind]) == expected_hash, f'source hash changed: {acc}/{label}/{kind}')
                check(fasta(root/source['fasta']) == genome, 'DNA/sequence IDs differ')
                models[label], raw, published = load_models(root, source)
                folder = output/'analyses'/(acc+'-'+label)
                inputs = ['-f', root/source['fasta'], '-g', root/source['gff']]
                run(acc+'-'+label+'-validate', ['validate', *inputs, '-o', output/'validation'/(acc+'-'+label), '--fail-on', 'error'])
                run(acc+'-'+label+'-analyze', ['analyze', *inputs, '-o', folder, '--table-format', 'tsv'])
                summaries[label] = check_analysis(folder, models[label], raw, genome)
                check_figures(folder, summary=summaries[label], models=models[label])
                check_analysis_terminal(output/'logs'/(acc+'-'+label+'-analyze.stdout'), summaries[label])
                analyses.append(folder/'summary.json')
                case_result['tools'][label] = dict(cds=len(models[label]), published_protein_checks=published,
                    computed_only_pseudogenes=len(models[label])-published, genome_bp=sum(map(len, genome.values())), sequence_count=len(genome))
            aggregate_args = ['compare']
            for label, source in sources.items():
                aggregate_args += ['--genome', label, root/source['fasta'], root/source['gff']]
            aggregate = output/'aggregate'/acc
            run(acc+'-aggregate', [*aggregate_args, '-o', aggregate])
            document = json_file(aggregate/'comparison.json')
            check('gene_level_comparison' not in document, 'gene comparison was not opt-in')
            check(len(document['datasets']) == len(sources) and len(document['pairwise_comparisons']) == math.comb(len(sources), 2), 'aggregate cardinality differs')
            check_figures(aggregate, profiles=document['datasets'])
            check_pairwise_metrics(document)
            for profile in document['datasets']:
                check_profile(profile, summaries[profile['label']], models[profile['label']])
            pair_jobs = []
            for left, right in itertools.combinations(sources, 2):
                base = output/'pairs'/acc/(left+'-vs-'+right)
                for a, b, direction in ((left, right, 'forward'), (right, left, 'reverse')):
                    command = ['compare', '--gene-level']
                    for label in (a, b):
                        command += ['--genome', label, root/sources[label]['fasta'], root/sources[label]['gff']]
                    pair_jobs.append((acc+'-'+a+'-'+b, [*command, '-o', base/direction]+(['-q'] if direction == 'reverse' else []), direction == 'reverse'))
            with ThreadPoolExecutor(max_workers=3) as executor:
                futures = [executor.submit(run, name, command, quiet) for name, command, quiet in pair_jobs]
                for future in futures:
                    future.result()
            for left, right in itertools.combinations(sources, 2):
                base = output/'pairs'/acc/(left+'-vs-'+right)
                forward = json_file(base/'forward/comparison.json')['gene_level_comparison']
                reverse = json_file(base/'reverse/comparison.json')['gene_level_comparison']
                terminal = {'Mode':'gene-level', 'Datasets':'2', 'Relationship rows':f"{len(forward['relationships']):,}",
                            'Accepted CDS pairs':f"{forward['status_counts'].get('matched',0):,}", 'Exact shared models':f"{forward['exact_model_count']:,}"}
                matched = forward['status_counts'].get('matched',0)
                for key, label in [('identical','Identical proteins'), ('terminal_only','Terminal differences'), ('initiator_recoding','Initiation context'), ('translation_difference','Other translation changes'), ('limited','Limited evidence')]:
                    count = forward['consequence_counts'].get(key,0)
                    terminal[label] = f'{count:,} ({100*count/matched:.2f}% of accepted pairs)' if matched else f'{count:,} (not available; no accepted pairs)'
                check_terminal(output/'logs'/(acc+'-'+left+'-'+right+'.stdout'), terminal)
                check(forward['cds_counts'] == {left:len(models[left]), right:len(models[right])}, 'gene per-source CDS counts differ')
                check_pair(forward, models[left], models[right])
                check_pair(reverse, models[right], models[left])
                groups = check_exports(base/'forward', forward)
                check_exports(base/'reverse', reverse)
                backward = {(r['id_b'], r['id_a']): r for r in reverse['relationships']}
                for row in forward['relationships']:
                    other = backward[row['id_a'], row['id_b']]
                    check(row['status'] == other['status'] and row['protein_consequence'] == other['protein_consequence'], 'swap status/consequence differs')
                    for key in ('start_delta_bp', 'stop_delta_bp', 'protein_length_delta_aa', 'n_terminal_delta_aa', 'c_terminal_delta_aa'):
                        check(other[key] == (-row[key] if row[key] is not None else None), 'swap signed delta differs: '+key)
                item = dict(a=left, b=right, rows=len(forward['relationships']), accepted=forward['status_counts'].get('matched', 0),
                            consequence_counts=forward['consequence_counts'], ambiguous_groups=groups,
                            report=str((base/'forward/comparison.html').relative_to(output)))
                case_result['pairs'].append(item)
                comparisons.append(item)
                print('AUDIT', acc, left, right, 'PASS', flush=True)
        except Exception as error:
            findings.append(dict(accession=acc, error=repr(error)))
            print('FAIL', acc, repr(error), flush=True)
        results_path = output/'verification.json'
        existing = json_file(results_path) if results_path.exists() else {'cases': []}
        existing['cases'] = [x for x in existing['cases'] if x['accession'] != acc]+[case_result]
        existing.update(findings=findings, limitation='Source consistency and computational verification; not experimental biological validation.')
        results_path.write_text(json.dumps(existing, indent=2)+'\n')
    if not findings:
        cohort = output/'cohort'
        run('cohort', ['summarize', *analyses, '-o', cohort])
        doc = json_file(cohort/'cohort.json')
        check_terminal(output/'logs/cohort.stdout', {'Samples':str(len(analyses))})
        check(doc['sample_count'] == len(analyses) and len(doc['samples']) == len(analyses), 'cohort lost/duplicated samples')
        with (cohort/'cohort.tsv').open(newline='') as f:
            table = list(csv.DictReader(f, delimiter='\t'))
        check(len(table) == len(doc['samples']), 'cohort TSV cardinality differs')
        for table_row, row in zip(table, doc['samples']):
            for key, value in table_row.items():
                expected = 'NA' if row[key] is None else f'{row[key]:.2f}' if isinstance(row[key], float) else str(row[key])
                check(value == expected, 'cohort TSV/JSON field differs: '+key)
        for row in doc['samples']:
            s = json_file(Path(row['summary_path']))
            check(row['cds_count'] == s['cds_count'] and row['scientific_fingerprint'] == s['scientific_fingerprint'], 'cohort join differs')
        # Three inputs must not silently produce a gene comparison.
        first = manifest['cases'][0]
        bad = ['compare', '--gene-level']
        for label, src in list(first['tools'].items())[:3]:
            bad += ['--genome', label, root/src['fasta'], root/src['gff']]
        run('reject-three-gene-inputs', [*bad, '-o', output/'rejected-three'], expected=2)
        check(not (output/'rejected-three/comparison.json').exists(), 'failed command published a report')
    if findings:
        raise SystemExit('FAIL: inspect '+str(output/'verification.json'))
    print(f'PASS {len(manifest["cases"])} assemblies, {len(analyses)} analyses, {len(comparisons)} distinct pairs and both directions; browser audit separate.', flush=True)


if __name__ == '__main__':
    main()
