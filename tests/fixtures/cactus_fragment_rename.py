"""Unchanged MIT-licensed helper from Cactus commit d7e5fc1b51327a5a4cdb494c85b0607e3d81a713.

Source: src/cactus/shared/common.py, get_faidx_subpath_rename_cmd.
Copyright Benedict Paten; upstream license: https://github.com/ComparativeGenomicsToolkit/cactus/blob/d7e5fc1b51327a5a4cdb494c85b0607e3d81a713/LICENSE.txt
Kept as a regression fixture, not imported by production code.
"""


def get_faidx_subpath_rename_cmd():
    return ['sed', '-e', r's/\([^:]*\):\([0-9]+\)-\([0-9]+\)/echo "\\1_sub_$((\\2-1))_\\3"/e']
