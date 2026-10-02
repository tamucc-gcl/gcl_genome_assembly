"""Compatibility entry point; production helper is the single implementation."""
from pathlib import Path
import runpy

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parents[2] / 'py_scripts' / 'sv_alignment_qc.py'), run_name='__main__')