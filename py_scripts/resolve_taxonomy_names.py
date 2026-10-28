"""Resolve exact scientific names in the run's local NCBI taxonomy snapshot."""
import argparse
import csv
import json
from pathlib import Path


def resolve(requests, names, nodes):
    normalize = lambda value: ' '.join(value.split()).casefold()
    queries = {normalize(row['species']) for row in requests}
    matches = {q: {} for q in queries}
    with Path(names).open(encoding='utf-8') as handle:
        for line in handle:
            fields = [v.strip() for v in line.split('|')]
            if len(fields) >= 4 and fields[3] == 'scientific name':
                key = normalize(fields[1])
                if key in matches:
                    matches[key][fields[0]] = fields[1]
    wanted = set().union(*matches.values())
    live = set()
    with Path(nodes).open(encoding='utf-8') as handle:
        for line in handle:
            taxid = line.split('|', 1)[0].strip()
            if taxid in wanted:
                live.add(taxid)
    result = []
    for row in requests:
        ids = matches[normalize(row['species'])]
        if not ids:
            raise ValueError(f"Unknown scientific name for {row['sample']}: {row['species']}; supply a valid taxid or scientific name")
        if len(ids) != 1:
            raise ValueError(f"Ambiguous scientific name for {row['sample']}: {row['species']}; supply a taxid")
        taxid = next(iter(ids))
        if taxid not in live:
            raise ValueError(f"Taxonomy snapshot has no node for matched taxid {taxid}")
        result.append(dict(sample=row['sample'], query=row['species'], species=ids[taxid], taxid=taxid))
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--requests', required=True)
    p.add_argument('--names', required=True)
    p.add_argument('--nodes', required=True)
    p.add_argument('--output', default='resolved_names.tsv')
    a = p.parse_args()
    try:
        result = resolve(json.loads(Path(a.requests).read_text()), a.names, a.nodes)
    except (ValueError, OSError) as exc:
        p.exit(65, str(exc) + '\n')
    with open(a.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=['sample', 'query', 'species', 'taxid'], delimiter='\t')
        writer.writeheader()
        writer.writerows(result)
