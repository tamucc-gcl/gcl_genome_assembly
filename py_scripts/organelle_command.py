"""Execute GetOrganelle without turning logged/zero-exit errors into no-recovery."""
import argparse
from pathlib import Path
import re
import subprocess
import sys


def validate_log(text, command=None):
    # GetOrganelle 1.7.7.1 catches exceptions/SystemExit and may return zero.
    if re.search(r'(?m)(?:^| - |\[)ERROR(?:\]|:)|^Traceback \(most recent call last\):',text):
        raise ValueError('GetOrganelle reported an error; this is not biological no-recovery')
    if 'Thank you!' not in text:
        # 1.7.7.1 from-assembly main() returns early (without its closing log)
        # when a successfully slimmed graph is empty. Verify that exact artifact
        # state; a zero exit or a partial log alone remains insufficient.
        if command and Path(command[0]).name == 'get_organelle_from_assembly.py' and '-o' in command:
            out = Path(command[command.index('-o') + 1])
            initial = out / 'initial_assembly_graph.fastg'
            slimmed = out / 'slimmed_assembly_graph.fastg'
            if (re.search(r'(?m) - INFO: Slimming .+ finished!\s*$', text)
                    and initial.is_file() and initial.stat().st_size > 0
                    and slimmed.is_file() and slimmed.stat().st_size == 0):
                for log in out.glob('*slim*.log*'):
                    if re.search(r'(?mi)ERROR|Traceback \(most recent call last\)',
                                 log.read_text(encoding='utf-8', errors='replace')):
                        raise ValueError('GetOrganelle slimming reported an error')
                return 'no_target_after_slimming'
        raise ValueError('GetOrganelle completion marker missing; run was not verified complete')
    return 'completed'


def run(command,log):
    with Path(log).open('wb') as output:
        process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        for block in iter(process.stdout.readline,b''):
            output.write(block);sys.stdout.buffer.write(block);sys.stdout.buffer.flush()
        code=process.wait()
    if code:raise SystemExit(code)
    outcome=validate_log(Path(log).read_text(encoding='utf-8',errors='replace'),command)
    if outcome == 'no_target_after_slimming':
        message='[organelle_command] Verified empty target graph after slimming; retry recovered no target contigs.\n'
        with Path(log).open('a',encoding='utf-8') as output:output.write(message)
        print(message,end='')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--log',required=True);p.add_argument('command',nargs=argparse.REMAINDER)
    a=p.parse_args();command=a.command[1:] if a.command[:1]==['--'] else a.command
    if not command:p.error('Missing GetOrganelle command')
    run(command,a.log)
