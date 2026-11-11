"""Independent checks for public annotation panels; no Annostat imports.

GFF3/output coordinates are 1-based inclusive; Biopython locations and internal
position arrays are 0-based. Source proteins establish annotation consistency,
not experimental correctness. Optional audit dependencies: gffutils, numpy.
"""
import csv
import json
import math
import re
from collections import Counter, defaultdict
from statistics import median

import gffutils
import numpy as np
from Bio import SeqIO
from Bio.Data import CodonTable

ATTRS = ('gene', 'product', 'Dbxref', 'pseudo', 'pseudogene', 'partial',
         'start_range', 'end_range', 'transl_except', 'exception')


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def fasta(path):
    result = {}
    for record in SeqIO.parse(path, 'fasta'):
        check(record.id not in result, 'duplicate FASTA ID: ' + record.id)
        result[record.id] = str(record.seq).upper()
    return result


def attr(attrs, key):
    if key not in attrs:
        return None
    text = ','.join(attrs[key])
    return ','.join(sorted(set(text.split(',')))) if key == 'Dbxref' else ' '.join(text.split())


def load_models(root, source):
    """Join independently parsed GFF, GenBank and published FAA by locus tag."""
    rows = list(gffutils.iterators.DataIterator(str(root / source['gff'])))
    parents = {attr(r.attributes, 'ID'): r for r in rows if r.featuretype in {'gene', 'pseudogene'}}
    circular = {r.seqid for r in rows if r.featuretype == 'region' and attr(r.attributes, 'Is_circular') == 'true'}
    groups = defaultdict(list)
    for r in rows:
        if r.featuretype == 'CDS':
            groups[attr(r.attributes, 'ID')].append(r)
    reference = {}
    genome = fasta(root/source['fasta'])
    seen_sequences = set()
    mapping = source.get('seqid_mapping', {})
    for record in SeqIO.parse(root / source['gbff'], 'genbank'):
        seqid = mapping.get(record.id, record.id)
        check(seqid in genome and str(record.seq).upper() == genome[seqid], 'GenBank/FASTA genome DNA differs')
        seen_sequences.add(seqid)
        for feature in record.features:
            if feature.type == 'CDS':
                key = (seqid, feature.qualifiers['locus_tag'][0])
                check(key not in reference, 'duplicate source locus tag')
                reference[key] = record, feature
    check(seen_sequences == set(genome), 'GenBank/FASTA sequence ID coverage differs')
    check(len(groups) == len(reference), 'GFF/GenBank CDS counts differ')
    published = fasta(root / source['faa'])
    result, used, published_ids = {}, set(), set()
    for identifier, pieces in groups.items():
        first = min(pieces, key=lambda r: r.start) if pieces[0].strand == '+' else max(pieces, key=lambda r: r.end)
        locus = attr(first.attributes, 'locus_tag') or identifier
        key = first.seqid, locus
        check(key in reference and key not in used, 'non-bijective GFF/GenBank join: ' + str(key))
        used.add(key)
        record, feature = reference[key]
        strand = '+' if feature.location.strand == 1 else '-'
        parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
        path = np.concatenate([np.arange(s, e, dtype=np.int32) if strand == '+' else np.arange(e-1, s-1, -1, dtype=np.int32) for s, e in parts])
        intervals = []
        for r in pieces:
            check(r.seqid == key[0] and r.strand == strand, 'GFF/GenBank strand differs')
            intervals.append((r.start-1, min(r.end, len(record))))
            if r.end > len(record):
                intervals.append((0, r.end-len(record)))
        check(sorted(parts) == sorted(intervals), 'GFF/GenBank intervals differ: ' + identifier)
        phase = int(feature.qualifiers.get('codon_start', ['1'])[0])-1
        check(phase == int(first.frame), 'GFF/GenBank initial phase differs')
        raw = dict(first.attributes)
        inherited = {}
        parent = parents.get(attr(raw, 'Parent'))
        if parent:
            inherited = {k: v for k, v in parent.attributes.items() if k in ATTRS}
            if parent.featuretype == 'pseudogene':
                inherited['pseudo'] = ['true']
        attrs = {**inherited, **raw}
        nt = str(feature.extract(record.seq)).upper()
        code = int(feature.qualifiers.get('transl_table', ['11'])[0])
        table = CodonTable.unambiguous_dna_by_id[code]
        marker = 'end_range' if strand == '-' else 'start_range'
        generic = attr(attrs, 'partial') in {'true', 'yes', '1'} and not ({'start_range', 'end_range'} & attrs.keys())
        complete = not phase and marker not in attrs and not generic
        coding = nt[phase:]
        if 'translation' in feature.qualifiers:
            protein = feature.qualifiers['translation'][0]
            protein_id = attr(raw, 'protein_id') or locus
            check(protein_id in published and published[protein_id].removesuffix('*') == protein,
                  'source FAA/GenBank protein mismatch: ' + protein_id)
            published_ids.add(protein_id)
        else:
            check('pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers,
                  'missing source translation without pseudogene context')
            aa = [table.forward_table.get(c, '*' if c in table.stop_codons else 'X')
                  for c in (coding[i:i+3] for i in range(0, len(coding)-2, 3))]
            if aa and complete and coding[:3] in table.start_codons:
                aa[0] = 'M'
            protein = ''.join(aa).removesuffix('*')
        positions = path[phase:phase+3*len(protein)]+1
        codons = [tuple(map(int, x)) for x in positions.reshape((-1, 3))]
        result[identifier] = dict(seqid=first.seqid, strand=strand, phase=phase,
            path=path, physical=np.unique(path), lo=int(path.min()), hi=int(path.max()),
            signature=(first.seqid, strand, phase, path.tobytes()),
            nt=nt, coding=coding, protein=protein, code=code, attrs=attrs, codons=codons,
            initiator=table.forward_table.get(coding[:3]) if complete and coding[:3] in table.start_codons else None,
            circular=first.seqid in circular, genome_length=len(record))
    check(used == set(reference), 'unused GenBank CDS records')
    check(published_ids == set(published), 'unused published FAA records')
    return result, rows, len(published_ids)


def relation(a, b):
    shared, i, j = np.intersect1d(a['path'], b['path'], return_indices=True)
    if a['signature'] == b['signature']:
        return 'exact', len(shared)
    frames = np.any((i-a['phase']) % 3 != (j-b['phase']) % 3)
    for m in (a, b):
        positions, counts = np.unique(m['path'], return_counts=True)
        for p in np.intersect1d(positions[counts > 1], shared):
            residues = set((np.flatnonzero(a['path'] == p)-a['phase']) % 3) | set((np.flatnonzero(b['path'] == p)-b['phase']) % 3)
            frames |= len(residues) != 1
    if frames:
        return 'frame_disagreement', len(shared)
    short, long = sorted((a['path'], b['path']), key=len)
    if a['path'][-1] == b['path'][-1] and np.array_equal(short, long[-len(short):]):
        return 'start_shift', len(shared)
    if a['path'][0] == b['path'][0] and np.array_equal(short, long[:len(short)]):
        return 'stop_shift', len(shared)
    return 'overlap', len(shared)


def expected_graph(a, b):
    """Exhaustive interval masks, independent of Annostat's sweep algorithm."""
    ids = list(b)
    low, high = np.array([b[k]['lo'] for k in ids]), np.array([b[k]['hi'] for k in ids])
    groups = np.array([(b[k]['seqid'], b[k]['strand']) for k in ids])
    links = {}
    for ai, am in a.items():
        mask = (low <= am['hi']) & (high >= am['lo']) & (groups[:, 0] == am['seqid']) & (groups[:, 1] == am['strand'])
        for j in np.flatnonzero(mask):
            kind, overlap = relation(am, b[ids[j]])
            if overlap:
                links[ai, ids[j]] = kind, overlap
    ca, cb = Counter(m['signature'] for m in a.values()), Counter(m['signature'] for m in b.values())
    common = ca.keys() & cb.keys()
    accepted = {(i, j) for (i, j), (kind, _) in links.items()
                if kind == 'exact' and ca[a[i]['signature']] == cb[b[j]['signature']] == 1}
    anchored_a, anchored_b = {i for i, _ in accepted}, {j for _, j in accepted}
    alternatives = {(i, j) for (i, j), (kind, _) in links.items() if kind in {'start_shift', 'stop_shift'}
                    and i not in anchored_a and j not in anchored_b and a[i]['signature'] not in common and b[j]['signature'] not in common}
    da, db = Counter(i for i, _ in alternatives), Counter(j for _, j in alternatives)
    accepted |= {(i, j) for i, j in alternatives if da[i] == db[j] == 1}
    anchored_a, anchored_b = {i for i, _ in accepted}, {j for _, j in accepted}
    remaining = {(i, j) for i, j in links if i not in anchored_a and j not in anchored_b}
    da, db = Counter(i for i, _ in remaining), Counter(j for _, j in remaining)
    expected = {pair: 'matched' for pair in accepted}
    expected.update({(i, j): 'candidate' if da[i] == db[j] == 1 else 'ambiguous' for i, j in remaining})
    expected.update({(i, None): 'unmatched' for i in a if i not in anchored_a and not da[i]})
    expected.update({(None, j): 'unmatched' for j in b if j not in anchored_b and not db[j]})
    return expected, links, sum(min(ca[s], cb[s]) for s in common)


def check_pair(gene, a, b):
    expected, links, exact = expected_graph(a, b)
    rows = gene['relationships']
    check(len(rows) == len(expected), 'relationship cardinality differs')
    check({(r['id_a'], r['id_b']): r['status'] for r in rows} == expected, 'relationship graph differs')
    check(gene['exact_model_count'] == exact, 'exact model count differs')
    check(gene['exact_agreement'] == 2*exact/(len(a)+len(b)), 'exact agreement denominator differs')
    for row in rows:
        ma, mb = a.get(row['id_a']), b.get(row['id_b'])
        for side, m in (('a', ma), ('b', mb)):
            if not m:
                continue
            expected_fields = {'cds_length': len(m['path']), 'protein_length': len(m['protein']),
                'phase': m['phase'], 'genetic_code': m['code'], 'start_codon': m['coding'][:3],
                'first_amino_acid': m['protein'][:1] or None,
                'first_codon_positions': ','.join(map(str, m['path'][m['phase']:m['phase']+3]+1))}
            for key, value in expected_fields.items():
                check(row[key+'_'+side] == value, f'{key} differs: {row["id_"+side]}')
            for key in ATTRS:
                check(row[key.lower()+'_'+side] == attr(m['attrs'], key), 'attribute differs: '+key)
        if ma and mb:
            kind, overlap = links[row['id_a'], row['id_b']]
            check((row['relationship'], row['overlap_bp']) == (kind, overlap), 'overlap/type differs')
            check(row['protein_identical'] == (ma['protein'] == mb['protein']), 'protein equality differs')
            direction = 1 if ma['strand'] == '+' else -1
            delta = [direction*int(mb['path'][k]-ma['path'][k]) for k in (0,-1)]
            if ma['circular'] or mb['circular']:
                if kind == 'exact':
                    delta = [0,0]
                elif kind in {'start_shift','stop_shift'}:
                    index = 0 if kind == 'start_shift' else 1
                    downstream = len(mb['path']) < len(ma['path']) if index == 0 else len(mb['path']) > len(ma['path'])
                    delta[index] = delta[index] % ma['genome_length'] if downstream else -((-delta[index]) % ma['genome_length'])
                    delta[1-index] = 0
                else:
                    delta = [None,None]
            check([row['start_delta_bp'],row['stop_delta_bp']] == delta, 'source-derived boundary offsets differ')
        if row['status'] != 'matched':
            check(row['protein_consequence'] == 'unresolved' and row['alignment_identity_percent'] is None,
                  'unaccepted link has a protein verdict')
            continue
        pa, pb, ka, kb = ma['protein'], mb['protein'], ma['codons'], mb['codons']
        if row['protein_consequence'] == 'limited':
            check(row['alignment_identity_percent'] is None, 'limited evidence has identity')
            check(set(pa+pb)-set('ACDEFGHIKLMNPQRSTVWYOU') or row['unapplied_transl_except_a'] or row['unapplied_transl_except_b'] or len(set(ka)) != len(ka) or len(set(kb)) != len(kb), 'unexplained limited protein evidence')
            continue
        check(row['protein_length_delta_aa'] == len(pb)-len(pa), 'protein delta differs')
        lookup = {key: j for j, key in enumerate(kb)}
        common = [(i, lookup[key]) for i, key in enumerate(ka) if key in lookup]
        check(common and all(y[1] > x[1] for x, y in zip(common, common[1:])), 'unordered common codons')
        same = initiation = other = 0
        for i, j in common:
            if pa[i] == pb[j]:
                same += 1
            elif (i == 0 and pa[i] == 'M' and pb[j] == ma['initiator']) or (j == 0 and pb[j] == 'M' and pa[i] == mb['initiator']):
                initiation += 1
            else:
                other += 1
        check((row['shared_codon_count'], row['identical_aa'], row['initiator_mismatch_count'], row['internal_mismatch_count']) == (len(common), same, initiation, other), 'shared-codon residue counts differ')
        left, right = (pa, pb) if row['protein_consequence'] == 'identical' else (row['alignment_a'], row['alignment_b'])
        check(left.replace('-', '') == pa and right.replace('-', '') == pb and len(left) == len(right), 'alignment does not retain source proteins')
        ia = ib = gaps = 0
        for x, y in zip(left, right):
            if x != '-' and y != '-':
                check(ka[ia] == kb[ib], 'aligned residues from different genomic codons')
            else:
                gaps += 1
            ia += x != '-'
            ib += y != '-'
        check(row['substitution_count'] == initiation+other and row['gap_aa'] == gaps, 'alignment counts differ')
        check(math.isclose(row['alignment_identity_percent'], 100*same/len(left)), 'identity denominator differs')
        fa, fb = common[0]
        la, lb = common[-1]
        check((row['n_terminal_delta_aa'], row['c_terminal_delta_aa']) == (fb-fa, (len(pb)-lb-1)-(len(pa)-la-1)), 'terminal deltas differ')
        for key, value in [('n_terminal_residues_a', pa[:fa]), ('n_terminal_residues_b', pb[:fb]), ('c_terminal_residues_a', pa[la+1:]), ('c_terminal_residues_b', pb[lb+1:])]:
            check(row[key] == value, 'terminal sequence differs')
        category = row['protein_consequence']
        if category == 'identical':
            check(pa == pb and ka == kb, 'false identical verdict')
        elif category == 'terminal_only':
            check(not (initiation or other) and gaps and row['relationship'] in {'start_shift', 'stop_shift'}, 'false terminal-only verdict')
        elif category == 'initiator_recoding':
            check(initiation and not other, 'false initiation-context verdict')
    labels = list(gene['cds_counts'])
    check(gene['unmatched_cds_counts'] == {label: sum(r['status'] == 'unmatched' and r['id_'+side] is not None for r in rows) for label, side in zip(labels, ('a','b'))}, 'unmatched CDS denominator differs')
    check(gene['status_counts'] == dict(Counter(r['status'] for r in rows)), 'status totals differ')
    check(gene['relationship_counts'] == dict(Counter(r['relationship'] for r in rows)), 'relationship totals differ')
    check(gene['consequence_counts'] == dict(Counter(r['protein_consequence'] for r in rows if r['status'] == 'matched')), 'consequence totals differ')


def check_analysis(folder, models, rows, genome):
    summary = json.loads((folder/'summary.json').read_text())
    nt, aa = fasta(folder/'sequences/cds_nucleotide.fasta'), fasta(folder/'sequences/cds_protein.fasta')
    check(set(nt) == set(aa) == set(models), 'exported CDS ID sets differ')
    for key, m in models.items():
        check(nt[key] == m['nt'] and aa[key] == m['protein'], 'exported sequence differs: '+key)
    length = sum(map(len, genome.values()))
    bases = Counter(''.join(genome.values()))
    gc = 100*(bases['G']+bases['C'])/sum(bases[x] for x in 'ACGT')
    covered = sum(len(np.unique(np.concatenate([m['physical'] for m in models.values() if m['seqid'] == sid]))) for sid in genome)
    check(summary['cds_count'] == len(models) and summary['genome_length'] == length, 'CDS/genome count differs')
    check(math.isclose(summary['quality_control']['genome_gc_percent'], gc), 'GC differs')
    check(math.isclose(summary['quality_control']['coding_density_percent'], 100*covered/length), 'union coding density differs')
    seen, types = set(), Counter()
    for i, r in enumerate(rows):
        key = (r.featuretype, attr(r.attributes, 'ID') or ('anonymous', i))
        if key not in seen:
            seen.add(key)
            types[r.featuretype] += 1
    check(summary['feature_type_counts'] == dict(types) and summary['gff_row_count'] == len(rows), 'feature row/model counts differ')
    check(summary['total_features'] == sum(types.values()), 'feature total differs')
    with (folder/'tables/features.tsv').open(newline='') as f:
        overview = list(csv.DictReader(f, delimiter='\t'))
    check([(r['ID'], r['Feature type'], int(r['Start position']), int(r['Stop position']), r['Strand']) for r in overview] ==
          [(attr(r.attributes,'ID') or '', r.featuretype, r.start, r.end, r.strand) for r in rows], 'feature TSV coordinates/IDs differ')
    check(summary['rna_counts'] == {k: v for k, v in types.items() if k.lower().endswith('rna')}, 'RNA count differs')
    categories, annotated = Counter(), 0
    for m in models.values():
        letters = [letter for token in m['attrs'].get('Dbxref', []) if token.startswith('COG:')
                   and token[4:].isalpha() and token[4:].isupper() for letter in token[4:]]
        categories.update(letters)
        annotated += bool(letters)
    check(summary['cog_category_counts'] == dict(categories) and summary['cds_with_cog_count'] == annotated, 'COG assignments/coverage differ')
    named = sum(bool(attr(m['attrs'], 'gene')) for m in models.values())
    hypothetical = sum('hypothetical protein' in (attr(m['attrs'], 'product') or '').lower() for m in models.values())
    check(summary['cds_with_gene_count'] == named and summary['hypothetical_cds_count'] == hypothetical, 'annotation coverage differs')
    starts = Counter(m['coding'][:3] for m in models.values() if len(m['coding']) >= 3)
    codons = Counter(m['coding'][i:i+3] for m in models.values() for i in range(0, len(m['coding'])-2, 3))
    check(summary['start_codon_counts'] == dict(starts) and summary['complete_codon_count'] == sum(codons.values()), 'codon totals differ')
    with (folder/'tables/codon_usage.csv').open() as f:
        table = list(csv.DictReader(f))
    check({r['codon']: int(r['count']) for r in table} == dict(codons), 'codon CSV differs')
    for r in table:
        check(abs(float(r['percentage'])-100*codons[r['codon']]/sum(codons.values())) < 1e-6, 'codon percentage differs')
    check(all((folder/x).is_file() for x in summary['output_files']), 'declared output missing')
    return summary


def check_figures(folder, summary=None, profiles=None, gene=None, models=None):
    """Read actual SVG bar geometry/labels, not the plotting functions."""
    from xml.etree import ElementTree as ET
    ns = '{http://www.w3.org/2000/svg}'
    def axes(path):
        return [g for g in ET.parse(path).iter(ns+'g') if g.get('id', '').startswith('axes_')]
    def texts(node):
        return [''.join(t.itertext()) for t in node.iter(ns+'text')]
    def scale(axis, letter):
        points = []
        for g in axis.iter(ns+'g'):
            if g.get('id', '').startswith(letter+'tick_'):
                t = next(g.iter(ns+'text'))
                try:
                    value = float(''.join(t.itertext()).replace(',', '').replace('−', '-'))
                    points.append((value, float(t.get(letter))))
                except (ValueError, TypeError):
                    continue
        check(len(points) >= 2, 'numeric SVG scale unavailable')
        a, b = points[0], points[-1]
        return (b[1]-a[1])/(b[0]-a[0]), a[1]-a[0]*(b[1]-a[1])/(b[0]-a[0])
    def bars(axis):
        result = []
        for g in axis:
            if g.tag != ns+'g' or not g.get('id', '').startswith('patch_'):
                continue
            for p in g.findall(ns+'path'):
                if p.get('clip-path') and 'fill: none' not in p.get('style', ''):
                    coords = list(map(float, re.findall(r'-?\d+(?:\.\d+)?(?:e[+-]?\d+)?', p.get('d', ''), re.I)))
                    if len(coords) == 8:
                        result.append(coords)
        return result
    if gene is not None:
        axis = axes(folder/'plots/gene_consequences.svg')[0]
        slope, _ = scale(axis, 'x')
        expected = [gene['consequence_counts'].get(k, 0) for k in ('identical', 'terminal_only', 'initiator_recoding', 'translation_difference', 'limited')]
        actual = [(v[2]-v[0])/slope for v in bars(axis)]
        check(np.allclose(actual, expected, atol=1e-3), 'gene SVG bar lengths differ')
        labels = Counter(texts(axis))
        required = Counter(f"{v:,} / {gene['status_counts'].get('matched', 0):,}" for v in expected)
        check(not (required-labels), 'gene SVG denominator labels differ')
        check('Accepted CDS pairs (count)' in labels, 'gene SVG terminology differs')
    if summary is not None:
        starts = summary['start_codon_counts']
        table = CodonTable.unambiguous_dna_by_id[summary['genetic_code']]
        expected = [starts.get(c, 0) for c in ('ATG', 'GTG', 'TTG')]
        expected += [sum(v for c, v in starts.items() if c in table.start_codons and c not in {'ATG','GTG','TTG'}),
                     sum(v for c, v in starts.items() if c not in table.start_codons), len(models)-sum(starts.values())]
        axis = axes(folder/'plots/start_codons.svg')[0]
        slope, _ = scale(axis, 'x')
        check(np.allclose([(v[2]-v[0])/slope for v in bars(axis)], expected, atol=1e-3), 'first-codon SVG bar lengths differ')
        for v in expected:
            check(f'{v:,} ({100*v/len(models):.2f}%)' in texts(axis), 'first-codon SVG labels differ')
        axis = axes(folder/'plots/cds_lengths.svg')[0]
        sx, ox = scale(axis, 'x')
        sy, _ = scale(axis, 'y')
        geometry = bars(axis)
        edges = [(v[0]-ox)/sx for v in geometry]+[(geometry[-1][2]-ox)/sx]
        lengths = np.array([len(m['path']) for m in models.values()])
        description = ET.parse(folder/'plots/cds_lengths.svg').getroot().find(ns+'desc').text
        trim = re.search(r'([\d,]+) CDS above ([\d,]+) nt', description)
        if trim:
            outliers, cutoff = (int(x.replace(',', '')) for x in trim.groups())
            check(cutoff == np.quantile(lengths, .99, method='lower') and np.sum(lengths > cutoff) == outliers, 'histogram cutoff disclosure differs')
            plotted = lengths[lengths <= cutoff]
        else:
            plotted = lengths
        counts, _ = np.histogram(plotted, bins=np.rint(edges).astype(int))
        heights = [abs(v[5]-v[1])/abs(sy) for v in geometry]
        check(np.allclose(heights, counts, atol=1e-3), 'CDS histogram heights differ')
        labels = texts(axis)
        check(f'Median: {median(lengths):,.0f} nt' in labels and f'Mean: {sum(lengths)/len(lengths):,.0f} nt' in labels, 'histogram summary markers differ')
        if summary['cog_data_available']:
            expected = [v for _, v in sorted(summary['cog_category_counts'].items(), key=lambda x: (-x[1], x[0]))]
            axis = axes(folder/'plots/cog_categories.svg')[0]
            slope, _ = scale(axis, 'x')
            check(np.allclose([(v[2]-v[0])/slope for v in bars(axis)], expected, atol=1e-3), 'COG SVG bar lengths differ')
    if profiles is not None:
        metrics = [('gc_percent','%'), ('coding_density_percent','%'), ('cog_coverage_percent','%'),
                   ('hypothetical_percent','%'), ('gene_name_percent','%'), ('cds_per_mb',' / Mb')]
        for axis, (key, suffix) in zip(axes(folder/'plots/dataset_overview.svg'), metrics):
            required = Counter('Not available' if p[key] is None else f'{p[key]:.2f}{suffix}' for p in profiles)
            check(not (required-Counter(texts(axis))), 'overview SVG measurement labels differ: '+key)
