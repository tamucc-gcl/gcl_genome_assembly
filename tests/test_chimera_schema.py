"""The empty/reference table must aggregate with populated/recovered tables."""
import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "py_scripts"))
from chimera_schema import CALLED_FIELDS, read_table, write_called
from break_chimeras import read_rows


class CalledJoinContractTests(unittest.TestCase):
    def test_empty_then_populated_aggregate_preserves_named_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            empty, populated = root / "empty.tsv", root / "populated.tsv"
            write_called(empty, [])
            write_called(populated, [dict(assembly="sample", scaffold="s",
                left_arm_fraction="0.4", structural_status="requires_adjudication",
                assembly_sha256="digest", coordinate_stage="pre_finishing")])
            self.assertEqual(read_table(empty)[0], read_table(populated)[0])
            combined = empty.read_text() + "".join(populated.read_text().splitlines(True)[1:])
            row = next(csv.DictReader(combined.splitlines(), delimiter="\t"))
            self.assertEqual(row["left_arm_fraction"], "0.4")
            self.assertEqual(row["assembly_sha256"], "digest")
            self.assertEqual(len(row), len(CALLED_FIELDS))

    def test_partial_review_row_is_normalized_without_shifting_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "recovered.tsv"
            write_called(path, [dict(assembly="a", gap_len="100", join_scope="recovered")])
            fields, rows = read_table(path)
            self.assertEqual(fields, CALLED_FIELDS)
            self.assertEqual(rows[0]["gap_len"], "100")
            self.assertEqual(rows[0]["left_arm_fraction"], ".")
            self.assertEqual(rows[0]["join_scope"], "recovered")

    def test_malformed_manual_instructions_fail_instead_of_disappearing(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "instructions.tsv"
            for text in ("assembly\tscaffold\na\ts\textra\n",
                         "assembly\tscaffold\na\n", "assembly\tassembly\na\tb\n"):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_table(path)
                with self.assertRaises(ValueError):
                    read_rows(str(path))


if __name__ == "__main__":
    unittest.main()
