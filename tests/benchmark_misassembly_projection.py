"""Reproduce the fragmented-label bottleneck; verify identical assignments."""
import argparse
import json
import time
from test_misassembly_performance import naive_tracks,compact_reference,tracks,blocks


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--blocks',type=int,default=3000)
    args=parser.parse_args();n=args.blocks
    if n<2:raise ValueError('At least two blocks required')
    length=11*n
    labels={'t':[dict(lo=11*i,hi=11*i+11,chrom='chr1' if i<n//2 else 'chr2') for i in range(n)]}
    paf=['q',str(length),'0',str(length),'+','t',str(length),'0',str(length),str(10*n),str(length),'60','cg:Z:'+'10=1X'*n]
    start=time.perf_counter();old=naive_tracks([paf],labels);old_seconds=time.perf_counter()-start
    start=time.perf_counter();new=tracks([paf],labels,mode='chromosome',min_mapq=20,minimum_target=0,compact=True);new_seconds=time.perf_counter()-start
    assert new==compact_reference(old),'Optimization changed assignments'
    print(json.dumps(dict(label_intervals=n,cigar_segments_before=len(blocks(paf)),cigar_segments_after=len(blocks(paf,compact=True)),
        output_intervals_before=len(old),output_intervals_after=len(new),previous_seconds=old_seconds,indexed_seconds=new_seconds,
        speedup=old_seconds/new_seconds,assignments_identical=True),indent=2))


if __name__=='__main__':main()
