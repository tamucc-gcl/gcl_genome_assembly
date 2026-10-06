import sys
import unittest
import tempfile
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'py_scripts'))
from chimera_adjudicate import decide
from chimera_sequence_context import narrow_spanning_molecules


class DecisionTests(unittest.TestCase):
    def complete(self):
        return dict(verified_gap=True, independent_discordant_individuals=3,
                    hifi_informative=True, hifi_spanning_molecules=0, hic_informative=True,
                    hic_support_loss=True, matched_controls_pass=True,
                    alternative_placements_checked=True, graph_contradiction=False)

    def test_every_required_assay_is_required(self):
        evidence = self.complete()
        self.assertEqual(decide(evidence)[:2], ('UNJOIN_UNSUPPORTED', True))
        for key in evidence:
            reduced = dict(evidence); del reduced[key]
            self.assertEqual(decide(reduced)[:2], ('UNRESOLVED', False), key)

    def test_noninformative_is_not_zero_support(self):
        for key in ('hifi_informative', 'hic_informative', 'matched_controls_pass', 'alternative_placements_checked'):
            evidence = self.complete(); evidence[key] = False
            self.assertFalse(decide(evidence)[1])

    def test_positive_continuity_veto_and_distinct_individual_threshold(self):
        self.assertEqual(decide({'hifi_spanning_molecules':2})[:2], ('RETAIN', False))
        evidence = self.complete(); evidence['independent_discordant_individuals'] = 2
        self.assertFalse(decide(evidence)[1])

    def test_deletion_bracketing_does_not_count_as_continuity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'reads.sam'
            path.write_text('good\t0\tscaffold\t1\t60\t4000M\t*\t0\t0\t*\t*\n'
                            'good\t0\tscaffold\t1\t60\t4000M\t*\t0\t0\t*\t*\n'
                            'deletion\t0\tscaffold\t1\t60\t1500M1000D1500M\t*\t0\t0\t*\t*\n'
                            'secondary\t256\tscaffold\t1\t60\t4000M\t*\t0\t0\t*\t*\n')
            self.assertEqual(narrow_spanning_molecules(path, 1500, 2500), 1)


if __name__ == '__main__':
    unittest.main()
