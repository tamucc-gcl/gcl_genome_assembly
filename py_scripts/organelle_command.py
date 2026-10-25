"""Execute GetOrganelle without turning logged/zero-exit errors into no-recovery."""
import argparse
from pathlib import Path
import re
import subprocess
import sys


def validate_log(text):
    # GetOrganelle 1.7.7.1 catches exceptions/SystemExit and may return zero.
    if re.search(r'(?m)(?:^| - |\[)ERROR(?:\]|:)|^Traceback \(most recent call last\):',text):
        raise ValueError('GetOrganelle reported an error; this is not biological no-recovery')
    if 'Thank you!' not in text:
        raise ValueError('GetOrganelle completion marker missing; run was not verified complete')


def run(command,log):
    with Path(log).open('wb') as output:
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        for block in iter(process.stdout.readline,b''):
            output.write(block);sys.stdout.buffer.write(block);sys.stdout.buffer.flush()
        code=process.wait()
    if code:raise SystemExit(code)
    validate_log(Path(log).read_text(encoding='utf-8',errors='replace'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log',required=True);p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:p.error('Missing GetOrganelle command')
    run(command,a.log)
