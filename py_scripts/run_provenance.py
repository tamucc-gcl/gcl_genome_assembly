"""Collect current task version evidence and resolved execution environments."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shlex


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_table(path,rows,fields):
    with Path(path).open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t');writer.writeheader();writer.writerows(rows)


def task_context(folder):
    folder=Path(folder)
    # Version outputs may be inside a task subdirectory (e.g. representation/).
    for parent in [folder]+list(folder.parents)[:3]:
        if (parent/'.command.run').is_file():folder=parent;break
    wrapper=folder/'.command.run';script=folder/'.command.sh'
    result=dict(workdir=str(folder),process='unavailable',wrapper_sha256=None,command_sha256=None,
                conda_prefixes=[],container_files=[],context_status='unavailable')
    if script.is_file():result['command_sha256']=sha256(script)
    if not wrapper.is_file():return result
    result['context_status']='recorded';result['wrapper_sha256']=sha256(wrapper)
    text=wrapper.read_text(errors='replace')
    task=re.search(r'^# NEXTFLOW TASK:\s*(.+)$',text,re.M)
    if task:result['process']=task.group(1)
    prefixes=set();images=set()
    for line in text.splitlines():
        try:tokens=shlex.split(line)
        except ValueError:continue
        for i,token in enumerate(tokens):
            if token.endswith('/bin/activate') and i+1<len(tokens):
                candidate=tokens[i+1]
                if Path(candidate).is_absolute():prefixes.add(candidate)
            if token.endswith(('.sif','.img','.simg')) and Path(token).is_absolute():images.add(token)
    result['conda_prefixes']=sorted(prefixes);result['container_files']=sorted(images)
    return result


def environment_inventory(contexts):
    conda={};containers={}
    for context in contexts:
        for prefix in context['conda_prefixes']:
            if prefix in conda:continue
            root=Path(prefix)/'conda-meta';records=[]
            if root.is_dir():
                for path in sorted(root.glob('*.json')):
                    data=json.loads(path.read_text())
                    records.append({key:data.get(key) for key in ['name','version','build','channel','url','sha256','md5']})
            conda[prefix]=dict(status='recorded' if records else 'unavailable',packages=records)
        for image in context['container_files']:
            if image in containers:continue
            path=Path(image)
            containers[image]=dict(status='recorded' if path.is_file() else 'unavailable',
                sha256=sha256(path) if path.is_file() else None,
                bytes=path.stat().st_size if path.is_file() else None)
    return dict(conda=conda,containers=containers)


def collect(specs,staged,context,out):
    out=Path(out);out.mkdir(exist_ok=True)
    files=sorted(Path(staged).glob('version*/*'))
    if len(files)!=len(specs):raise ValueError('Version source/staging cardinality mismatch')
    rows=[];issues=[];sources=[];contexts=[]
    # Numbered staging directories preserve the input-list order; lexical sorting
    # alone would put version100 before version11.
    files.sort(key=lambda p:int(re.search(r'(\d+)$',p.parent.name).group(1)))
    for spec,path in zip(specs,files):
        ctx=task_context(Path(spec['path']).parent);contexts.append(ctx)
        source=dict(path=spec['path'],sha256=sha256(path),task=ctx);sources.append(source)
        lines=path.read_text(errors='replace').splitlines()
        observed=0
        for number,line in enumerate(lines,1):
            if not line.strip():continue
            fields=line.split('\t')
            if fields in [['tool','version'],['process','tool','version']]:continue
            if len(fields)==2:tool,version=fields;process=ctx['process']
            elif len(fields)==3:process,tool,version=fields
            else:
                issues.append(dict(source=spec['path'],line=number,issue='unlabeled_or_malformed_record',value=line));continue
            if not tool.strip() or not version.strip():
                issues.append(dict(source=spec['path'],line=number,issue='empty_tool_or_version',value=line));continue
            rows.append(dict(tool=tool.strip(),version=version.strip(),process=process,source=spec['path']))
            observed+=1
        if not observed:issues.append(dict(source=spec['path'],line=0,issue='no_labeled_version_records',value=''))
    rows=sorted({tuple(r.values()) for r in rows})
    write_table(out/'software_versions.tsv',[dict(zip(['tool','version','process','source'],r)) for r in rows],
                ['tool','version','process','source'])
    write_table(out/'software_version_issues.tsv',issues,['source','line','issue','value'])
    payload=dict(context=context,sources=sources,version_records=len(rows),issues=len(issues),
        scope='Version-emitting tasks; cached producers retain their original task context. Not a claim of complete task coverage.')
    (out/'run_provenance.json').write_text(json.dumps(payload,indent=2)+'\n')
    (out/'runtime_environments.json').write_text(json.dumps(environment_inventory(contexts),indent=2)+'\n')
    catalog=[];definitions=sorted(Path('catalog').glob('source*/*'),key=lambda p:int(re.search(r'(\d+)$',p.parent.name).group(1)))
    if len(definitions)!=len(context['catalog_sources']):raise ValueError('Dependency catalog staging mismatch')
    for original,path in zip(context['catalog_sources'],definitions):
        digest=sha256(path);catalog.append(dict(source=original,sha256=digest,kind='source_file',line=0,specification='',pin_policy=''))
        in_dependencies=False
        for number,line in enumerate(path.read_text().splitlines(),1):
            stripped=line.strip()
            if stripped=='dependencies:':in_dependencies=True;continue
            if path.suffix in ['.yaml','.yml'] and in_dependencies and stripped.startswith('- '):
                spec=stripped[2:].split(' #')[0].strip()
                policy='version_constraint' if any(c in spec for c in '=<>') else 'unconstrained'
                catalog.append(dict(source=original,sha256=digest,kind='conda_dependency',line=number,specification=spec,pin_policy=policy))
            elif path.suffix in ['.nf','.config'] and not stripped.startswith('//') and re.search(r'\b(?:\w+_)?(?:container|conda)\s*=?\s*[\'\"]',line):
                policy='mutable_reference' if re.search(r'(:master|:latest|/latest/)',line) else 'inspect_declaration'
                catalog.append(dict(source=original,sha256=digest,kind='execution_declaration',line=number,specification=stripped,pin_policy=policy))
    write_table(out/'dependency_catalog.tsv',catalog,['source','sha256','kind','line','specification','pin_policy'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['specs','staged','context','out']:p.add_argument('--'+name,required=True)
    a=p.parse_args();collect(json.loads(Path(a.specs).read_text()),a.staged,json.loads(Path(a.context).read_text()),a.out)


if __name__=='__main__':main()
