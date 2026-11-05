#!/usr/bin/env python3
"""Validate published report files after Nextflow has finished publishing."""
import argparse
import csv
import html
import json
import os
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


def rows(path):
    with path.open(encoding='utf-8') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def anchors(text):
    result=set(re.findall(r'<a\b[^>]*(?:id|name)=["\']([^"\']+)',text,re.I))
    counts={}
    for heading in re.findall(r'^#{1,6}\s+(.+?)\s*#*$',text,re.M):
        heading=re.sub(r'<[^>]+>','',html.unescape(heading)).lower()
        key=re.sub(r'[^\w -]','',heading).replace(' ','-')
        number=counts.get(key,0);counts[key]=number+1
        result.add(key+('-'+str(number) if number else ''))
    return result


def validate(root):
    root=Path(os.path.abspath(root))
    errors=[];checked=0
    manifest=root/'reports/main/report_manifest.json'
    # A run with reporting disabled is explicitly not a report-validation success.
    provenance=root/'provenance/run_provenance.json'
    if provenance.is_file():
        controls=json.loads(provenance.read_text(encoding='utf-8')).get('context',{}).get('parameters',{})
        if controls.get('run_report') is False:
            return dict(status='DISABLED',errors=[],checked_links=0)
    if not manifest.is_file():
        return dict(status='FAIL',errors=['Main report manifest missing'],checked_links=0)
    checks=root/'reports/main/report_completeness.tsv'
    if not checks.is_file():errors.append('Report completeness checks missing')
    else:
        errors.extend(r['Entity']+': expected '+r['Feature']+' missing' for r in rows(checks) if r['Status']=='MISSING')
    registry=root/'reports/main/artifact_registry.tsv'
    artifact_docs=[]
    if not registry.is_file():errors.append('Artifact registry missing')
    else:
        for row in rows(registry):
            path=root/row['path']
            if not path.is_file() or path.stat().st_size==0:errors.append('Registered artifact missing or empty: '+row['path'])
            if path.suffix.lower()=='.md':artifact_docs.append(path)
    docs=list(dict.fromkeys([root/'assembly_report.md',*(root/'reports/main').glob('*.md'),*artifact_docs]))
    queued=set(docs)
    if not docs[0].is_file():errors.append('Main Markdown report missing')
    for doc in docs:
        if not doc.is_file():continue
        text=doc.read_text(encoding='utf-8')
        # Generated reports use inline Markdown links and HTML image tags.
        links=[(value,marker=='!') for marker,value in re.findall(r'(!?)\[[^\]\n]*\]\(([^)\n]+)\)',text)]
        links += [(value,attribute.lower()=='src') for attribute,value in re.findall(r'\b(src|href)=["\']([^"\']+)["\']',text,re.I)]
        for value,is_image in set(links):
            value=html.unescape(value.strip().strip('<>'))
            url=urlsplit(value)
            if url.scheme or url.netloc:continue
            # Published symlinks intentionally point into Nextflow work directories.
            # Validate their published URL location without resolving the symlink.
            target=Path(os.path.abspath(doc.parent/unquote(url.path))) if url.path else doc
            if not target.is_relative_to(root):
                errors.append(str(doc.relative_to(root))+': link escapes results directory: '+value);continue
            checked+=1
            if not is_image and target.is_dir():
                continue  # Directory navigation is a valid report link.
            if not target.is_file() or target.stat().st_size==0:
                errors.append(str(doc.relative_to(root))+': missing/empty target '+value);continue
            if target.suffix.lower()=='.md' and target not in queued:
                queued.add(target);docs.append(target)
            if url.fragment and target.suffix.lower()=='.md' and unquote(url.fragment) not in anchors(target.read_text(encoding='utf-8')):
                errors.append(str(doc.relative_to(root))+': unknown heading '+value)
    return dict(status='FAIL' if errors else 'PASS',errors=sorted(set(errors)),checked_links=checked)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();result=validate(args.root)
    if args.root.is_dir():
        # Keep post-run diagnostics outside published task directories, which may
        # be symlinks into the Nextflow cache and must not be modified afterward.
        (args.root/'report_validation.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (args.root/'report_validation.md').write_text('# Published report validation\n\n'+result['status']+'\n\n'+str(result['checked_links'])+' local links checked. This validates report outputs and links, not biological acceptance.\n\n'+'\n'.join('- '+e for e in result['errors'])+'\n',encoding='utf-8')
    print(result['status']+': published report validation ('+str(result['checked_links'])+' links)')
    for error in result['errors']:print(error)
    raise SystemExit(1 if result['status']=='FAIL' else 0)


if __name__=='__main__':main()
