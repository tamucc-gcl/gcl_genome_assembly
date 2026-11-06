"""Bounded native fixtures: paired assemblies, complete HiFi reads and haploid skip."""
import json
from pathlib import Path
import random
import sys


def setup(directory):
    import pysam
    root=Path(directory);root.mkdir(parents=True,exist_ok=True)
    rng=random.Random(125)
    ref=''.join(rng.choices('ACGT',k=250000))
    alternate=list(ref)
    for i in range(150,249900,173):alternate[i]=rng.choice([b for b in 'ACGT' if b!=ref[i]])
    alternate=''.join(alternate)
    own_read_reference=list(ref)
    own_read_reference[2000]=next(b for b in 'ACGT' if b!=own_read_reference[2000])
    for name,sequence in [('hap1',ref),('hap2',alternate),('read_reference',''.join(own_read_reference))]:
        p=root/(name+'.fa')
        p.write_text('>long_scaffold\n'+sequence+'\n>short_piece\n'+'A'*1000+'\n')
        pysam.faidx(str(p))
    fq=root/'complete.fastq'
    with fq.open('w') as f:
        for i in range(1500):
            start=rng.randrange(len(ref)-5000);seq=(ref if rng.random()<.5 else alternate)[start:start+5000]
            f.write(f'@read_{i}\n{seq}\n+\n'+('I'*len(seq))+'\n')
    rows=[]
    for sample,ploidy,haps in [('assembly_ind',2,['hap1','hap2']),('read_ind',2,['primary']),('haploid_ind',1,['primary'])]:
        for hap in haps:
            source='read_reference' if sample=='read_ind' else 'hap2' if hap=='hap2' else 'hap1'
            rows.append(dict(id=sample+'_'+hap,sample=sample,taxid='1',species='Fixture species',ploidy=ploidy,
                             haplotype=hap,representation='hic_phased' if len(haps)==2 else 'collapsed_or_primary',
                             completion='finalized',fasta=str(root/(source+'.fa')),fai=str(root/(source+'.fa.fai'))))
    (root/'eligibility.json').write_text(json.dumps(dict(assemblies=rows)))
    (root/'manifest.tsv').write_text('sample\tmode\nassembly_ind\tassemblies\nread_ind\thifi\nhaploid_ind\tauto\n')
    (root/'reads.json').write_text(json.dumps([dict(sample='read_ind',fastq=str(fq))]))
    return root


if __name__=='__main__':setup(sys.argv[1])
