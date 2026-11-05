"""Render actual retained junctions without changing source evidence or mappings."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'py_scripts'))
from junction_plot_data import export
from misassembly_core import table


def main():
    p=argparse.ArgumentParser();p.add_argument('--assemblies',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--rscript',required=True)
    p.add_argument('--actions',type=Path)
    a=p.parse_args();script=Path(__file__).resolve().parents[1]/'r_scripts/plot_junction_placements.R'
    n=0
    for packet in sorted(a.assemblies.iterdir()):
        if not (packet/'evidence.json').is_file():continue
        data=json.loads((packet/'evidence.json').read_text())
        if not any(e.get('junction_evidence',{}).get('hifi_status')=='measured' for e in data['events']):continue
        target=a.out/packet.name;target.mkdir(parents=True,exist_ok=True)
        shutil.copy2(packet/'evidence.json',target/'evidence.json')
        shutil.copytree(packet/'junction/hifi',target/'junction/hifi',dirs_exist_ok=True)
        for e in data['events']:
            if export(target,e):
                cuts=[r for r in table(a.actions) if r['assembly']==data['assembly'] and r['scaffold']==e['scaffold']
                      and r['event_id']==e['id'] and r['selected']=='YES'] if a.actions else []
                selected=target/'selected_cuts.json';selected.write_text(json.dumps(cuts))
                subprocess.run([a.rscript,str(script),str(target),e['id'],str(selected)],check=True)
                for ext in ('png','pdf'):
                    assert (target/'junction'/(e['id']+'.placements.'+ext)).stat().st_size>1000
                n+=1
    assert n>0,'No actual junction assays rendered'
    print(f'Rendered {n} retained junctions; source evidence unchanged')


if __name__=='__main__':main()
