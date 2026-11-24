"""A fake runtime rejects any attempt to run an image shell and records app arguments."""
from pathlib import Path


def setup(root):
    root=Path(root)
    tools=root/'tools'
    tools.mkdir(parents=True,exist_ok=True)
    binary=tools/'singularity'
    binary.write_text('''#!/usr/bin/env python3
import json,sys
from pathlib import Path
a=sys.argv[1:]
if a==['--version']:
    print('singularity fixture 1');raise SystemExit
if '/bin/bash' in a or '/bin/sh' in a:
    print('Shell requested inside shell-less image',file=sys.stderr);raise SystemExit(91)
assert a[0]=='exec',a
app=a.index('/app/fcs/bin/av_screen_x')
assert a[app+1:app+4]==['-o','/output-volume','--euk'],a
assert '--cleanenv' in a and '--pwd' in a,a
binds=[a[i+1] for i,x in enumerate(a) if x=='--bind']
output=Path(next(x[:-len(':/output-volume')] for x in binds if x.endswith(':/output-volume')))
source=Path(next(x[:-len(':/sample-volume:ro')] for x in binds if x.endswith(':/sample-volume:ro')))
assert (source/Path(a[-1]).name).is_file(),a
output.mkdir(exist_ok=True)
(output/'fcs_adaptor_report.txt').write_text('#accession\\tlength\\taction\\trange\\tname\\n')
(output/'runtime_arguments.json').write_text(json.dumps(a))
''')
    binary.chmod(0o755)
    return binary


if __name__=='__main__':setup('/tmp/fcs-adaptor-host-check')
