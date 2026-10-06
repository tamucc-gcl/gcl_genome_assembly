#!/usr/bin/env python3
"""Inventory retained hifiasm context without mapping, graph reconstruction or edits."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil


def trace_task_roots(trace, work_root, sample):
    """Resolve only matching hifiasm hashes, never recurse over the work tree."""
    roots, audit = [], []
    with trace.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'name', 'hash'}.issubset(reader.fieldnames or []):
            raise ValueError('Trace requires name and hash columns')
        for task in reader:
            name = task['name']
            if not name.endswith(f':HIFIASM ({sample})'):
                continue
            match = re.fullmatch(r'([0-9a-f]{2})/([0-9a-f]+)', task['hash'])
            if not match:
                raise ValueError('Invalid hifiasm task hash')
            parent = work_root/match[1]
            matches = sorted(p for p in parent.glob(match[2]+'*') if p.is_dir())
            status = 'resolved' if len(matches) == 1 else 'missing' if not matches else 'ambiguous'
            audit.append(dict(task=name, hash=task['hash'], status=status,
                              matches=[str(p.resolve()) for p in matches]))
            if len(matches) == 1 and matches[0] not in roots:
                roots.append(matches[0])
    return roots, audit


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--task-dir', type=Path, action='append', default=[],
                   help='Optional known hifiasm task directory; no recursive work-directory search')
    p.add_argument('--trace', type=Path, help='Resolve this sample hifiasm task from retained pipeline trace')
    p.add_argument('--work-root', type=Path, help='Work root used by that trace; required with --trace')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', a.sample) or a.sample in ('.', '..'):
        p.error('Sample must be a plain identifier')
    if bool(a.trace) != bool(a.work_root):
        p.error('--trace and --work-root must be supplied together')
    a.out.mkdir(parents=True, exist_ok=False)
    roots = [a.results/'assembly/contig/hifiasm'] + a.task_dir
    task_audit = []
    if a.trace:
        traced, task_audit = trace_task_roots(a.trace, a.work_root, a.sample)
        roots.extend(p for p in traced if p not in roots)
    records = []
    for i, root in enumerate(roots):
        files = sorted(root.glob(a.sample+'*')) if root.is_dir() else []
        files += [root/name for name in ('versions.tsv', '.command.sh', '.command.log')
                  if (root/name).is_file()]
        for path in files:
            if not path.is_file():
                continue
            stat = path.stat()
            record = dict(root_index=i, path=str(path.resolve()), name=path.name,
                          bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
            # Keep small diagnostic text; do not copy or hash gigabytes of graphs,
            # FASTA or caches in a context-presence check.
            if path.name.endswith(('.log', '.bed')) or path.name in ('versions.tsv', '.command.sh'):
                if stat.st_size <= 2000000:
                    dest = a.out/f'root_{i}'; dest.mkdir(exist_ok=True)
                    shutil.copy2(path, dest/path.name)
                    record['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            records.append(record)
    names = {r['name'] for r in records}
    graph_names = sorted(n for n in names if n.endswith('.gfa'))
    status = dict(sample=a.sample, roots=[str(r.resolve()) for r in roots], files=records,
                  task_resolution=task_audit,
                  graphs=graph_names,
                  raw_unitig_graph_present=any(n.endswith('.r_utg.gfa') for n in names),
                  processed_unitig_graph_present=any(n.endswith('.p_utg.gfa') for n in names),
                  hifi_overlap_caches_present=all(any(n.endswith(s) for n in names)
                      for s in ('.ec.bin', '.ovlp.source.bin', '.ovlp.reverse.bin')),
                  limitation='Presence and sizes only; no sequence/cache compatibility or path correctness validation.',
                  cut_authorized=False)
    (a.out/'context_inventory.json').write_text(json.dumps(status, indent=2)+'\n')
    print(f"Inventoried {len(records)} files; raw/processed unitig graphs: "
          f"{status['raw_unitig_graph_present']}/{status['processed_unitig_graph_present']}", flush=True)


if __name__ == '__main__':
    main()
