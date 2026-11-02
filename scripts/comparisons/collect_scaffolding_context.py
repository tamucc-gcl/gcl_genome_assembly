#!/usr/bin/env python3
"""Read-only task-context collection for comparing two scaffolding runs.

Diagnostic only: does not run Nextflow, modify its cache, or copy reads/assemblies.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import tarfile


def collect(old_trace, new_trace, assembly, work, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    records = []
    with tarfile.open(output, 'w:gz') as archive:
        def put(name, data):
            if isinstance(data, str):
                data = data.encode()
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))

        for label, trace in [('old', old_trace), ('new', new_trace)]:
            put(f'{label}/trace.tsv', trace.read_bytes())
            with trace.open() as handle:
                tasks = list(csv.DictReader(handle, delimiter='\t'))
            for ordinal, task in enumerate(tasks):
                name = task['name']
                selected = (assembly in name and name.startswith(
                    ('HIC_SCAFFOLDING:', 'CONTIG_REFINEMENT:')))
                if not selected and 'FIND_MITO_REFERENCE' not in name:
                    continue
                entry = {'run': label, 'task': task, 'files': [], 'missing': []}
                directory = Path(task.get('workdir', '-'))
                if not directory.is_dir():
                    prefix = task['hash'].split('/')
                    matches = list((work / prefix[0]).glob(prefix[1] + '*'))
                    directory = matches[0] if len(matches) == 1 else None
                entry['resolved_workdir'] = str(directory) if directory else None
                records.append(entry)
                if directory is None:
                    entry['missing'].append('task directory')
                    continue
                base = f'{label}/{ordinal:04d}'
                candidates = [directory / n for n in
                              ['.command.sh', '.command.out', '.command.err', 'versions.tsv']]
                candidates += list(directory.glob('*.agp'))
                if (directory / 'agp').is_dir():
                    candidates += list((directory / 'agp').glob('*.agp'))
                for path in dict.fromkeys(candidates):
                    if not path.is_file():
                        entry['missing'].append(str(path.relative_to(directory)))
                        continue
                    # Explicitly retain log tails, rather than copying giant tool logs.
                    limit = 2_000_000
                    size = path.stat().st_size
                    if size > limit and path.suffix == '.agp':
                        entry['missing'].append(f'{path.name}: over {limit} bytes')
                        continue
                    with path.open('rb') as handle:
                        if size > limit:
                            handle.seek(size-limit)
                        data = handle.read(limit)
                    member = base + '/' + str(path.relative_to(directory)).replace('\\', '/')
                    put(member, data)
                    entry['files'].append({'path': member, 'source_bytes': size,
                                           'tail_only': size > limit})
                links = {}
                for folder in [directory, directory/'assembly', directory/'agp', directory/'pairs']:
                    if folder.is_dir():
                        for path in folder.iterdir():
                            if path.is_symlink():
                                links[str(path.relative_to(directory))] = str(path.readlink())
                put(base+'/staging_links.json', json.dumps(links, indent=2))
        put('context_inventory.json', json.dumps(records, indent=2))
    print(f'Collected {len(records)} task records: {output}')
    missing = sum('task directory' in r['missing'] for r in records)
    if missing:
        print(f'WARNING: {missing} task directories unavailable; listed in inventory')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('old_trace', type=Path)
    parser.add_argument('new_trace', type=Path)
    parser.add_argument('--assembly', required=True)
    parser.add_argument('--work', type=Path, default=Path('work'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    collect(args.old_trace, args.new_trace, args.assembly, args.work, args.output)
