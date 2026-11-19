"""Synthetic CLI/report audit for description-summary scope and interpretation.

Uses only tiny constructed fixtures; coordinates are 1-based inclusive.
Run: python examples/check_description_summary.py --output tmp/description-cases
"""
import argparse
import csv
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import quote

from Bio.Seq import Seq


def cds(start, end, strand='+', phase=0, extra='', identifier='{id}'):
    return f'chr\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\tID={identifier};product={{product}};gene={{gene}}{extra}'


def cases():
    rows = []
    def add(name, dna, a, b=None, category='identical', proteins=None, accepted=1,
            unresolved=(0, 0), products=('DUF1234 protein', 'hypothetical protein'), labels=('A', 'B')):
        rows.append(dict(name=name, dna=dna, a=a, b=a if b is None else b,
                         category=category, proteins=proteins, accepted=accepted,
                         unresolved=unresolved, products=products, labels=labels))
    add('complete-plus', 'ATGAAATAA', [cds(1,9)], proteins=('MK','MK'))
    add('complete-minus', str(Seq('ATGAAATAA').reverse_complement()), [cds(1,9,'-')], proteins=('MK','MK'))
    for phase in (1,2):
        dna='A'*phase+'ATGAAATAA'
        add('phase-'+str(phase), dna, [cds(1,len(dna),phase=phase)], proteins=('MK','MK'))
    add('multipart-plus', 'AATGA'+'C'*15+'AATAA',
        [cds(1,5,phase=1),cds(21,25,phase=2)], proteins=('MK','MK'))
    add('multipart-minus', 'TTATT'+'A'*15+'TCATT',
        [cds(21,25,'-',1),cds(1,5,'-',2)], proteins=('MK','MK'))
    region='chr\ttest\tregion\t1\t18\t.\t+\t.\tID=region;Is_circular=true'
    circular='AAATAA'+'C'*9+'ATG'
    add('circular-plus',circular,[region,cds(16,24)],proteins=('MK','MK'))
    add('circular-minus',str(Seq(circular).reverse_complement()),[region,cds(13,21,'-')],proteins=('MK','MK'))
    add('table-4','ATGTGACAATAA',[cds(1,12,extra=';transl_table=4')],proteins=('MWQ','MWQ'))
    add('table-25','ATGTGACAATAA',[cds(1,12,extra=';transl_table=25')],proteins=('MGQ','MGQ'))
    add('different-tables','ATGTGACAATAA',[cds(1,12,extra=';transl_table=4')],
        [cds(1,12,extra=';transl_table=25')],category='translation_difference',proteins=('MWQ','MGQ'))
    add('partial-both','GTGAAATAA',[cds(1,9,extra=';partial=true')],proteins=('VK','VK'))
    add('partial-one','GTGAAATAA',[cds(1,9)],[cds(1,9,extra=';partial=true')],
        category='initiator_recoding',proteins=('MK','VK'))
    for amino, code in [('Sec','U'),('Pyl','O')]:
        add('supported-'+amino,'ATGTGACAATAA',[cds(1,12,extra=f';transl_except=(pos:4..6,aa:{amino})')],
            proteins=('M'+code+'Q','M'+code+'Q'))
    add('unapplied-exception','ATGAAAGGGTAA',[cds(1,12,extra=';transl_except=(pos:join(4..5,6),aa:Sec)')],
        category='limited',proteins=('MKG','MKG'))
    add('unknown-residue','ATGNNNAAATAA',[cds(1,12)],category='limited',proteins=('MXK','MXK'))
    add('internal-stop','ATGTGACAATAA',[cds(1,12)],category='limited',proteins=('M*Q','M*Q'))
    add('same-first-codon-different-position','ATGGCTATGGAACAACAATAA',[cds(1,21)],[cds(7,21)],
        category='terminal_only',proteins=('MAMEQQ','MEQQ'))
    parent='chr\ttest\tgene\t1\t9\t.\t+\t.\tID={id}-parent;product={product};gene={gene}'
    child='chr\ttest\tCDS\t1\t9\t.\t+\t0\tID={id};Parent={id}-parent'
    add('parent-attributes','ATGAAATAA',[parent,child],proteins=('MK','MK'))
    add('missing-products','ATGAAATAA',[cds(1,9)],proteins=('MK','MK'),products=(None,'  '))
    add('uncertain-other-text','ATGAAATAA',[cds(1,9)],proteins=('MK','MK'),
        products=('uncharacterized UPF1234 protein','unknown protein'))
    add('escaped-metadata','ATGAAATAA',[cds(1,9)],proteins=('MK','MK'),
        products=('</script><img src=x onerror=alert(1)> μ','hypothetical protein'),
        labels=('A <& μ>', 'B "quoted"'))
    add('duplicate-exact','ATGAAATAA',[cds(1,9,identifier='{id}1'),cds(1,9,identifier='{id}2')],
        [cds(1,9)],accepted=0,unresolved=(2,1))
    add('ambiguous-overlaps','ATG'*40,[cds(1,30,identifier='{id}1'),cds(31,60,identifier='{id}2')],
        [cds(16,45,identifier='{id}1'),cds(46,75,identifier='{id}2')],accepted=0,unresolved=(2,2))
    add('frame-candidate','ATG'*10,[cds(1,30)],[cds(1,30,phase=1)],accepted=0,unresolved=(1,1))
    add('opposite-strands','ATG'*10,[cds(1,30)],[cds(1,30,'-')],accepted=0,unresolved=(1,1))
    add('one-empty-source','ATGAAATAA',[cds(1,9)],[],accepted=0,unresolved=(1,0))
    add('no-cds','ATGAAATAA',[],[],accepted=0)
    return rows


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def run(command, log, quiet=False, expected=0):
    result = subprocess.run([sys.executable,'-m','annostat',*map(str,command)],capture_output=True,text=True,timeout=120)
    Path(str(log)+'.stdout').write_text(result.stdout,encoding='utf-8')
    Path(str(log)+'.stderr').write_text(result.stderr,encoding='utf-8')
    require(result.returncode == expected, f'{log.name}: exit {result.returncode}, expected {expected}: {result.stderr[-1000:]}')
    require('Traceback' not in result.stderr+result.stdout, 'unexpected traceback')
    if quiet and expected == 0:
        require(not result.stdout and not result.stderr,'quiet mode wrote routine output')
    if expected:
        require(result.stderr.strip() and not result.stdout.strip(),'invalid input stream contract')
    return result


def check_result(folder, spec, reverse):
    saved=json.loads((folder/'comparison.json').read_text(encoding='utf-8'))
    gene=saved['gene_level_comparison']
    description=gene['description_summary']
    labels=list(spec['labels'][::-1] if reverse else spec['labels'])
    accepted=spec['accepted']
    identical=accepted if spec['category']=='identical' else 0
    require(gene['status_counts'].get('matched',0)==accepted,'accepted fixture count')
    require(gene['consequence_counts']==({spec['category']:accepted} if accepted else {}),'fixture protein category')
    require(description['accepted_pairs']==accepted and description['identical_protein_pairs']==identical,'fixture summary denominator')
    unresolved=spec['unresolved'][::-1] if reverse else spec['unresolved']
    require(description['unresolved_cds_counts']==dict(zip(labels,unresolved)),'distinct unresolved fixture models')
    if accepted:
        row=next(r for r in gene['relationships'] if r['status']=='matched')
        if spec['proteins']:
            proteins=spec['proteins'][::-1] if reverse else spec['proteins']
            for side, protein in zip(('a','b'),proteins):
                require(row['protein_length_'+side]==len(protein),'fixture predicted protein length')
                require(row['first_amino_acid_'+side]==protein[:1],'fixture first residue')
                if row['alignment_'+side] is not None:
                    require(row['alignment_'+side].replace('-','')==protein,'fixture retained protein sequence')
    cells=description['product_status_matrix']
    require(len(cells)==9 and sum(c['pair_count'] for c in cells)==identical,'fixture matrix partition')
    if identical:
        expected_status={'missing-products':('missing','missing'),
                         'uncertain-other-text':('other_product_text','other_product_text')}.get(
                             spec['name'],('other_product_text','explicit_hypothetical'))
        if reverse: expected_status=expected_status[::-1]
        nonzero=[c for c in cells if c['pair_count']]
        require(len(nonzero)==1 and (nonzero[0]['product_status_a'],nonzero[0]['product_status_b'])==expected_status,'fixture matrix cell')
    else:
        require(all(c['percent_of_identical_pairs'] is None for c in cells),'empty denominator became measured zero')
    with (folder/'tables/gene_product_status.tsv').open(encoding='utf-8',newline='') as handle:
        table=list(csv.DictReader(handle,delimiter='\t'))
    require(len(table)==9,'fixture TSV matrix size')
    for cell, exported in zip(cells,table):
        require(exported['pair_count']==str(cell['pair_count']),'fixture TSV count')
        value=cell['percent_of_identical_pairs']
        require(exported['percent_of_identical_pairs']==('' if value is None else str(value)),'fixture TSV missingness')
    html=(folder/'comparison.html').read_text(encoding='utf-8')
    payload=json.loads(re.search(r'<script id="gene-data" type="application/json">(.*?)</script>',html,re.S).group(1))
    require(payload==gene,'fixture HTML/JSON payload')
    require('</script><img' not in html,'unescaped source markup')
    require('1-based inclusive' in html,'fixture coordinate disclosure')
    return {'case':spec['name'],'direction':'reverse' if reverse else 'forward','accepted':accepted,
            'identical':identical,'category':spec['category'] if accepted else 'unresolved'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=True)
    results=[]
    for spec in cases():
        base=root/spec['name']
        base.mkdir(exist_ok=True)
        fasta=base/'genome.fna'
        fasta.write_text('>chr\n'+spec['dna']+'\n',encoding='utf-8')
        paths=[]
        for index, side in enumerate(('a','b')):
            gff=base/(side+'.gff3')
            product=spec['products'][index]
            text='##gff-version 3\n'+'\n'.join(
                row.replace('{id}',side).replace('{product}',quote(product or '',safe='')).replace('{gene}','symbol_'+side)
                for row in spec[side])+'\n'
            gff.write_text(text,encoding='utf-8')
            paths.append((spec['labels'][index],fasta,gff))
        for reverse in (False,True):
            output=base/('reverse' if reverse else 'forward')
            command=['compare','--gene-level']
            for label,fna,gff in paths[::-1] if reverse else paths:
                command+=['--genome',label,fna,gff]
            command+=['-o',output]
            if reverse: command+=['--quiet']
            result=run(command,base/('reverse' if reverse else 'forward'),quiet=reverse)
            results.append(check_result(output,spec,reverse))
            if not reverse:
                require(re.search(r'^\s*Accepted CDS pairs\s+'+str(spec['accepted'])+r'\s*$',result.stdout,re.M),'fixture terminal accepted count')
        print('PASS',spec['name'],flush=True)
        (root/'verification.json').write_text(json.dumps({'cases':results},indent=2)+'\n',encoding='utf-8')
    # A label-only perturbation must not affect matching or protein measurements.
    original=json.loads((root/'complete-plus/forward/comparison.json').read_text(encoding='utf-8'))['gene_level_comparison']
    renamed=json.loads((root/'escaped-metadata/forward/comparison.json').read_text(encoding='utf-8'))['gene_level_comparison']
    for key in ('status_counts','relationship_counts','consequence_counts','exact_agreement'):
        require(original[key]==renamed[key],'label-only control changed '+key)
    (root/'verification.json').write_text(json.dumps({'cases':results,'label_only_control':'passed'},indent=2)+'\n',encoding='utf-8')
    print(f'PASS {len(cases())} synthetic examples in both directions',flush=True)


if __name__=='__main__':
    main()
