#!/usr/bin/env python3
"""Inventory retained hifiasm context without mapping, graph reconstruction or edits."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--task-dir', type=Path, action='append', default=[],
                   help='Optional known hifiasm task directory; no recursive work-directory search')
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', a.sample) or a.sample in ('.', '..'):
        p.error('Sample must be a plain identifier')
    a.out.mkdir(parents=True, exist_ok=False)
    roots = [a.results/'assembly/contig/hifiasm'] + a.task_dir
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
