import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts/comparisons'))
from prepare_shortread_smoke import prepare, record


class ShortreadSmokeTests(unittest.TestCase):
    def test_malformed_fastq_rejected(self):
        with self.assertRaises(ValueError):
            record(io.StringIO('@r/1\nAC\n+\n!\n'))

    def test_pair_prefix_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for mate in (1, 2):
                (root / str(mate)).write_text(f'@a/{mate}\nAC\n+\n!!\n@b/{mate}\nGT\n+\n!!\n')
            out = root / 'out'
            prepare(root / '1', root / '2', out, 1, 'smoke', 1)
            with gzip.open(out / 'reads_1.fastq.gz', 'rt') as handle:
                self.assertEqual(handle.read(), '@a/1\nAC\n+\n!!\n')
            self.assertEqual(json.loads((out / 'subset.json').read_text())['written_pairs'], 1)
            with self.assertRaises(FileExistsError):
                prepare(root / '1', root / '2', out, 1, 'smoke', 1)

    def test_mismatched_pair_has_no_completion_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / '1').write_text('@a/1\nAC\n+\n!!\n')
            (root / '2').write_text('@b/2\nAC\n+\n!!\n')
            with self.assertRaises(ValueError):
                prepare(root / '1', root / '2', root / 'out', 1, 'smoke', 1)
            self.assertFalse((root / 'out/subset.json').exists())


if __name__ == '__main__':
    unittest.main()
