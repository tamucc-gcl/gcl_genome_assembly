#!/usr/bin/env python3
"""Run NCBI's shell-less adaptor image directly from a host Nextflow task."""
import argparse
from collections import deque
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess


def launch(fasta, output, image, mode, engine, log, provenance):
    fasta = Path(fasta).resolve(strict=True)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if mode not in ('euk', 'prok'):
        raise ValueError('FCS-adaptor mode must be euk or prok')
    if any(c in str(fasta.parent) + str(output) for c in (':', ',')):
        raise ValueError('Container bind paths must not contain colons or commas')
    if engine in ('singularity', 'apptainer'):
        runtime = shutil.which(engine)
        if not runtime and engine == 'singularity':
            runtime = shutil.which('apptainer')
        if not runtime:
            raise RuntimeError('Singularity/Apptainer is required on the host running the adaptor task')
        if '://' not in image:
            image = str(Path(image).resolve(strict=True))
        command = [runtime, 'exec', '--cleanenv', '--bind', f'{fasta.parent}:/sample-volume:ro',
                   '--bind', f'{output}:/output-volume', '--pwd', '/output-volume', image]
    elif engine == 'docker':
        runtime = shutil.which('docker')
        if not runtime:
            raise RuntimeError('Docker is required on the host running the adaptor task')
        image = image.removeprefix('docker://')
        if '://' in image or Path(image).is_file():
            raise ValueError('Docker requires a registry image reference; use Singularity for a local SIF')
        command = [runtime, 'run', '--rm', '--init', '--user', f'{os.getuid()}:{os.getgid()}',
                   '--volume', f'{fasta.parent}:/sample-volume:ro', '--volume', f'{output}:/output-volume',
                   '--workdir', '/output-volume', image]
    else:
        raise ValueError('Container engine must be singularity, apptainer or docker')
    command += ['/app/fcs/bin/av_screen_x', '-o', '/output-volume', '--'+mode,
                '/sample-volume/'+fasta.name]
    version = subprocess.check_output([runtime, '--version'], text=True).strip()
    with Path(log).open('w') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        with Path(log).open() as handle:
            tail = ''.join(deque(handle, maxlen=40))
        raise RuntimeError(f'FCS-adaptor exited {result.returncode}; log: {log}\n{tail}')
    if not (output/'fcs_adaptor_report.txt').is_file():
        raise RuntimeError('FCS-adaptor did not create fcs_adaptor_report.txt; see '+str(log))
    digest = hashlib.sha256()
    with fasta.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8*1024*1024), b''):
            digest.update(chunk)
    Path(provenance).write_text(json.dumps(dict(image=image, runtime=version, mode=mode,
        source_sha256=digest.hexdigest(), command=command,
        launch='host runtime executes /app/fcs/bin/av_screen_x directly; no shell inside image'), indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('fasta','output','log','provenance'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--image',required=True)
    parser.add_argument('--mode',choices=['euk','prok'],required=True)
    parser.add_argument('--engine',choices=['singularity','apptainer','docker'],required=True)
    args=parser.parse_args()
    launch(args.fasta,args.output,args.image,args.mode,args.engine,args.log,args.provenance)
