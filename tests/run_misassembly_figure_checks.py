"""Check rendered figures against retained evidence, including missing-data routes."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import re
from urllib.parse import unquote,urlsplit
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from misassembly_core import table,sha
from misassembly_plot import detail_window


def check_actual(source,rendered,rscript=None):
    count=0
    for path in source.glob('*/evidence.json'):
        data=json.loads(path.read_text())
        target=rendered/path.parent.name
        if not target.is_dir():continue
        for e in data['events']:
            if rscript:
                before=sha(target/'evidence.json')
                script=Path(__file__).resolve().parents[1]/'r_scripts/plot_misassembly_evidence.R'
                subprocess.run([rscript,str(script),str(target),e['id'],str(target/'plots'/(e['id']+'.selected_cuts.json'))],check=True)
                assert sha(target/'evidence.json')==before,'Renderer modified evidence'
            audit={r['track']:int(r['rows']) for r in table(target/'plots'/(e['id']+'.render_audit.tsv'))}
            expected=dict(assigned_peer_intervals=sum(bool(b.get('chrom')) for t in e['peer_tracks'] for b in t['blocks']),
                consensus_intervals=len(e.get('chromosome_regions',[])),
                depth20=len(e.get('depth',[])),depth_all_primary=len(e.get('depth_all_primary',[])),
                contact_profile=sum(len(x) for x in e.get('contact_profiles',{}).values()),telomeres=len(e.get('telomeres',[])))
            for key,name in [('hifi20','hifi'),('hifi30','hifi_strict'),('hifi_all_primary','hifi_all_primary')]:
                expected[key]=len(e.get(name+'_profile',e.get(name,{})).get('probes',[]))
            assert audit==expected,(e['id'],audit,expected)
            selected=json.loads((target/'plots'/(e['id']+'.selected_cuts.json')).read_text())
            spec=json.loads((target/'plots'/(e['id']+'.render_spec.json')).read_text())
            assert spec['detail_window_bp']==list(detail_window(e,selected))
            for ending in ('.png','.pdf','.detail.png','.detail.pdf'):
                assert (target/'plots'/(e['id']+ending)).stat().st_size>1000
            count+=1
    assert count>0
    linked=0
    for doc in rendered.parent.rglob('*.md'):
        for url in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',doc.read_text(encoding='utf-8')):
            target=urlsplit(url)
            if target.scheme or target.netloc or not target.path:continue
            path=doc.parent/unquote(target.path)
            assert path.exists(),(doc,url)
            if path.is_file():assert path.stat().st_size>0,(doc,url)
            linked+=1
    print(f'PASS: {count} actual events retain every track interval/profile row and accepted detail window')
    print(f'PASS: {linked} local links/images in detailed reports resolve')


def synthetic(rscript):
    script=Path(__file__).resolve().parents[1]/'r_scripts/plot_misassembly_evidence.R'
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);(root/'plots').mkdir();selected=root/'selected.json'
        selected.write_text('[{"cut_bp":"700000"}]')
        e=dict(id='T1',scaffold='s',lo=500000,hi=600000,left='chr14',right='chr4',peer_tracks=[],cut_options=[],
            hifi=dict(probes=[dict(position=500000,molecules=10),dict(position=600000,molecules=8)]),
            hifi_all_primary=dict(probes=[dict(position=500000,molecules=400),dict(position=600000,molecules=8)]))
        (root/'evidence.json').write_text(json.dumps(dict(assembly='synthetic',events=[e])))
        subprocess.run([rscript,str(script),str(root),'T1',str(selected)],check=True)
        spec=json.loads((root/'plots/T1.render_spec.json').read_text())
        assert spec['detail_window_bp']==[600000,800000]
        assert spec['spanning_scale']['limit']<20
        assert spec['focal_colours']=={'chr4':'#0072B2','chr14':'#D55E00'}
        audit={r['track']:int(r['rows']) for r in table(root/'plots/T1.render_audit.tsv')}
        assert audit['contact_profile']==audit['depth20']==audit['telomeres']==0
        assert (root/'plots/T1.png').stat().st_size>1000
    print('PASS: missing assays stay unavailable; primary scales resist ambiguous spikes; user-cut close-up and focal colours preserved')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rscript',required=True)
    p.add_argument('--source',type=Path);p.add_argument('--rendered',type=Path);p.add_argument('--render',action='store_true');a=p.parse_args()
    synthetic(a.rscript)
    if a.source:check_actual(a.source,a.rendered,a.rscript if a.render else None)
