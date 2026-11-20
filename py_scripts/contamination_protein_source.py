#!/usr/bin/env python3
"""Resolve a pinned UniProt protein source and retain verified downloads for rebuilding."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

BASE = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/'


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def catalog_request(text, release):
    document = ET.fromstring(text)
    ns = {'m': 'http://www.metalinker.org/'}
    observed = document.findtext('m:version', namespaces=ns)
    if observed != release:
        raise ValueError(f'UniProt catalog release is {observed}, requested {release}; '
                         'supply an archived source and SHA256 for an older release')
    name = f'Reference_Proteomes_{release}.tar.gz'
    matches = [entry for entry in document.findall('m:files/m:file', ns) if entry.get('name') == name]
    if len(matches) != 1:
        raise ValueError('Expected one canonical reference-proteome archive in the UniProt catalog')
    entry = matches[0]
    md5 = entry.findtext('m:verification/m:hash[@type="md5"]', namespaces=ns)
    size = entry.findtext('m:size', namespaces=ns)
    if not re.fullmatch(r'[a-fA-F0-9]{32}', md5 or '') or not (size or '').isdigit() or int(size) < 1:
        raise ValueError('UniProt catalog lacks a usable archive checksum/size')
    return dict(release=release, source=BASE + name, source_format='tar', md5=md5.lower(),
                sha256=None, bytes=int(size), catalog=BASE + 'RELEASE.metalink')


def resolve(root, release, source=None, sha256=None, source_format='tar', catalog_text=None):
    if not re.fullmatch(r'[a-zA-Z0-9][a-zA-Z0-9_.-]*', release):
        raise ValueError('Invalid protein database release label')
    if bool(source) != bool(sha256):
        raise ValueError('A custom protein source requires both its source and SHA256')
    if sha256 and not re.fullmatch(r'[a-fA-F0-9]{64}', sha256):
        raise ValueError('Invalid expected protein source SHA256')
    if source and not source.startswith('https://'):
        source = str(Path(source).resolve())
    bundle = root / 'diamond/uniprot_reference_proteomes' / release
    database, manifest = bundle / 'reference_proteomes.dmnd', bundle / 'protein_database.json'
    if database.is_file() and manifest.is_file():
        identity = json.loads(manifest.read_text())
        if (identity.get('release') != release or identity.get('collection') != 'uniprot_reference_proteomes'
                or identity.get('source_format') != source_format or database.stat().st_size == 0
                or not re.fullmatch(r'[a-f0-9]{64}', identity.get('source_sha256', ''))):
            raise ValueError('Stored protein database identity differs from the requested release/format')
        if source:
            if identity.get('source') != source or identity['source_sha256'] != sha256.lower():
                raise ValueError('Stored protein database differs from the requested custom source')
        elif identity.get('source') != BASE + f'Reference_Proteomes_{release}.tar.gz':
            raise ValueError('Stored database is not the requested official UniProt archive')
        return dict(existing=True, bundle=str(bundle), release=release)
    if source:
        if not source.startswith('https://') and not Path(source).is_file():
            raise ValueError('Protein source file does not exist: ' + source)
        request = dict(release=release, source=source, sha256=sha256.lower(), md5=None,
                       source_format=source_format, bytes=None)
    else:
        if source_format != 'tar' or not re.fullmatch(r'\d{4}_\d{2}', release):
            raise ValueError('Official UniProt mode requires a YYYY_NN release and tar format')
        if catalog_text is None:
            catalog_text = subprocess.check_output(
                ['curl', '--fail', '--location', '--retry', '3', BASE + 'RELEASE.metalink'], text=True)
        request = catalog_request(catalog_text, release)
    return dict(request, existing=False, bundle=str(bundle))


def file_digests(path):
    sha256, md5 = hashlib.sha256(), hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            sha256.update(block)
            md5.update(block)
    return dict(sha256=sha256.hexdigest(), md5=md5.hexdigest())


def fetch(request, output):
    output.mkdir(parents=True, exist_ok=True)
    data = output / 'archive'
    if request['source'].startswith('https://'):
        partial = output / 'archive.partial'
        size = request.get('bytes')
        if size:
            remaining = max(0, size - (partial.stat().st_size if partial.exists() else 0))
            reserve = 10 * 1024**3 if size > 1024**3 else 10 * 1024**2
            if shutil.disk_usage(output).free < remaining + reserve:
                raise ValueError(f'Insufficient shared storage for {remaining} download bytes plus reserve')
        # Preserve partial bytes between attempts; verify the entire file before publication.
        for attempt in range(4):
            result = subprocess.run(['curl', '--fail', '--location', '--connect-timeout', '60',
                                     '--speed-time', '600', '--speed-limit', '1024', '--continue-at', '-',
                                     '--output', str(partial), request['source']])
            if result.returncode == 0:
                break
        else:
            raise RuntimeError('Protein archive download failed after four resumable attempts')
        actual = file_digests(partial)
        size = partial.stat().st_size
    else:
        partial = Path(request['source'])
        actual = file_digests(partial)
        size = partial.stat().st_size
    for name in ('sha256', 'md5'):
        if request.get(name) and actual[name] != request[name]:
            raise ValueError('Protein source ' + name.upper() + ' mismatch')
    if request.get('bytes') is not None and size != request['bytes']:
        raise ValueError('Protein archive size differs from the pinned catalog')
    if partial.parent == output:
        partial.rename(data)
    else:
        # A local, user-owned source is reused in place, with no extra archive-sized copy.
        data.symlink_to(partial.resolve())
    write_json(output / 'source.json', dict(request, **{k: v for k, v in actual.items()},
                                           actual_bytes=size))


def validate_download(request, data, manifest):
    identity = json.loads(manifest.read_text())
    for field in ('release', 'source', 'source_format'):
        if request[field] != identity.get(field):
            raise ValueError('Stored protein download differs from requested ' + field)
    for field in ('sha256', 'md5'):
        if request.get(field) and request[field] != identity.get(field):
            raise ValueError('Stored protein download checksum differs from the request')
    if data.stat().st_size != identity.get('actual_bytes'):
        raise ValueError('Stored protein archive size changed')
    return identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('resolve')
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--release', required=True)
    p.add_argument('--source')
    p.add_argument('--sha256')
    p.add_argument('--format', choices=['fasta', 'tar'], default='tar')
    p.add_argument('--out', type=Path, required=True)
    f = commands.add_parser('fetch')
    f.add_argument('--request', type=Path, required=True)
    f.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'resolve':
        write_json(args.out, resolve(args.root.resolve(), args.release, args.source, args.sha256, args.format))
    else:
        fetch(json.loads(args.request.read_text()), args.out)


if __name__ == '__main__':
    main()
