"""Verify one individual's selected final assembly for PSMC, independently of graphs."""
import argparse
import json
from pathlib import Path
from pangenome_input_audit import inspect_fasta


def verify(selection, fasta, fai):
    digest, lengths = inspect_fasta(fasta)
    indexed = {}
    for line in Path(fai).read_text(encoding='utf-8').splitlines():
        fields = line.split('\t')
        if len(fields) < 2 or fields[0] in indexed or int(fields[1]) <= 0:
            raise ValueError('Invalid individual PSMC reference FAI')
        indexed[fields[0]] = int(fields[1])
    if lengths != indexed:
        raise ValueError('Individual PSMC reference FASTA/FAI disagree')
    return dict(selection, sha256=digest, coordinate_stage='final',
        sequences=[dict(name=n, length=s) for n, s in lengths.items()],
        bases=sum(lengths.values()), sequence_count=len(lengths),
        ranking_basis='Deterministic hap1 for paired assemblies; own sole assembly otherwise. No QC superiority inferred.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['selection','fasta','fai']:p.add_argument('--'+name,required=True)
    a=p.parse_args()
    result=verify(json.loads(Path(a.selection).read_text(encoding='utf-8')),a.fasta,a.fai)
    Path('reference.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    Path('reference.md').write_text('# Individual PSMC reference\n\n'
        f"Individual: **{result['sample']}**. Reference: **{result['id']}**.\n\n"
        f"Selection: `{result['selection']}`. SHA256: `{result['sha256']}`.\n\n"
        'The final FASTA and FAI agree. Selection is deterministic within the individual; '
        'it does not establish biological correctness or rank haplotypes by quality.\n\n'
        '[Reference identity and sequence dictionary](reference.json)\n',encoding='utf-8')


if __name__=='__main__':main()
