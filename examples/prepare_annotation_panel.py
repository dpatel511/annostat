"""Prepare the checksum-pinned public annotation panel (test data only)."""
import argparse
import datetime
import gzip
import hashlib
import json
import shutil
import socket
import subprocess
import tarfile
import urllib.request
from pathlib import Path

from Bio import SeqIO
import pyrodigal


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprints(path):
    result = {}
    for record in SeqIO.parse(path, 'fasta'):
        if record.id in result:
            raise ValueError('duplicate FASTA ID: '+record.id)
        result[record.id] = [len(record), hashlib.sha256(bytes(record.seq).upper()).hexdigest()]
    return result


def download(info, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        temporary = Path(str(path)+'.part')
        if 'www.homd.org' in info['url']:
            # Try advertised IPv4 addresses, preserving TLS hostname verification.
            addresses = sorted({r[4][0] for r in socket.getaddrinfo('www.homd.org', 443, type=socket.SOCK_STREAM)})
            for address in addresses:
                result = subprocess.run(['curl', '--resolve', 'www.homd.org:443:'+address,
                    '--fail', '--location', '--connect-timeout', '10', '--max-time', '90', '-sS',
                    info['url'], '-o', str(temporary)], capture_output=True, text=True)
                if result.returncode == 0:
                    break
            else:
                raise RuntimeError('HOMD download failed: '+result.stderr)
        else:
            with urllib.request.urlopen(info['url'], timeout=90) as source, temporary.open('wb') as out:
                shutil.copyfileobj(source, out)
        temporary.replace(path)
    if path.stat().st_size != info['bytes'] or sha(path) != info['sha256']:
        raise ValueError('Frozen source mismatch: '+str(path)+'; do not replace expected hashes automatically')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path = Path(__file__).with_name('annotation-panel-manifest.json')
    manifest = json.loads(manifest_path.read_text())
    root = args.output.resolve()
    if pyrodigal.__version__ != '3.7.1':
        parser.error('this frozen panel requires test-only pyrodigal==3.7.1')
    for case in manifest['cases']:
        for info in case['downloads']:
            path = root/info['path']
            download(info, path)
            if path.suffix == '.gz' and not path.name.endswith('.tar.gz'):
                with gzip.open(path, 'rb') as source, path.with_suffix('').open('wb') as out:
                    shutil.copyfileobj(source, out)
        acc, tools = case['accession'], case['tools']
        archive = next(root/x['path'] for x in case['downloads'] if x['path'].endswith('.tar.gz'))
        with tarfile.open(archive) as package:
            for extension, key in [('fna', 'fasta'), ('gff3', 'gff'), ('gbff', 'gbff'), ('faa', 'faa')]:
                members = [m for m in package.getmembers() if m.isfile() and Path(m.name).name == acc+'.'+extension]
                if len(members) != 1:
                    raise ValueError('archive member not unique: '+acc+'.'+extension)
                path = root/tools['Bakta'][key]
                path.parent.mkdir(parents=True, exist_ok=True)
                with package.extractfile(members[0]) as source, path.open('wb') as out:
                    shutil.copyfileobj(source, out)
        canonical = fingerprints(root/tools['RefSeq']['fasta'])
        if canonical != case['sequences'] or fingerprints(root/tools['Bakta']['fasta']) != canonical:
            raise ValueError('assembly sequence mismatch: '+acc)
        if 'Prokka' in tools:
            source = tools['Prokka']
            mapping = source['seqid_mapping']
            original = fingerprints(root/source['source_fasta'])
            if len(original) != len(canonical) or {mapping[k]: v for k, v in original.items()} != canonical:
                raise ValueError('identifier mapping requires exact DNA identity: '+acc)
            with (root/source['source_gff']).open() as src, (root/source['gff']).open('w') as out:
                for line in src:
                    if line.strip() == acc.replace('GCF_', 'GCA_')+'|##gff-version 3':
                        line = '##gff-version 3\n'
                    elif line.startswith('##sequence-region '):
                        fields = line.split(); fields[1] = mapping[fields[1]]
                        line = ' '.join(fields)+'\n'
                    elif line.startswith('##FASTA'):
                        raise ValueError('unexpected embedded FASTA')
                    elif line.strip() and not line.startswith('#'):
                        fields = line.rstrip('\n').split('\t')
                        if len(fields) != 9:
                            raise ValueError('invalid GFF row')
                        fields[0] = mapping[fields[0]]
                        line = '\t'.join(fields)+'\n'
                    out.write(line)
        source = tools['Pyrodigal']
        finder = pyrodigal.GeneFinder(meta=False, closed=False)
        records = list(SeqIO.parse(root/source['fasta'], 'fasta'))  # Confirmed small panel: <=13.1 Mb per assembly.
        finder.train(*(bytes(r.seq) for r in records), translation_table=11)
        for key in ('gff', 'faa', 'gbff'):
            (root/source[key]).parent.mkdir(parents=True, exist_ok=True)
        with (root/source['gff']).open('w') as gff, (root/source['faa']).open('w') as aa, (root/source['gbff']).open('w') as gb:
            for i, record in enumerate(records):
                genes = finder.find_genes(bytes(record.seq))
                genes.write_gff(gff, record.id, header=i == 0, include_translation_table=True, full_id=True)
                genes.write_translations(aa, record.id, full_id=True)
                genes.write_genbank(gb, record.id, date=datetime.date.fromisoformat(manifest['prediction_date']))
        for label, source in tools.items():
            for kind, digest in source['sha256'].items():
                if sha(root/source[kind]) != digest:
                    raise ValueError(f'prepared hash mismatch: {acc}/{label}/{kind}')
        print('PREPARED', acc, ', '.join(tools), flush=True)
    (root/'input-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
