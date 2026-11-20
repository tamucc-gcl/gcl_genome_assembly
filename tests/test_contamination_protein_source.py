import hashlib
import gzip
import io
import json
from pathlib import Path
import sys
import tempfile
import tarfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'py_scripts'))
from contamination_protein_source import BASE, catalog_request, fetch, resolve, validate_download
from build_contamination_protein_db import prepare


def catalog(release='2026_04', digest='1' * 32):
    return f'''<metalink xmlns="http://www.metalinker.org/"><version>{release}</version><files>
    <file name="Reference_Proteomes_{release}.tar.gz"><size>1234</size>
    <verification><hash type="md5">{digest}</hash></verification></file></files></metalink>'''


class ProteinSourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_official_catalog_pins_release_size_and_checksum(self):
        request = catalog_request(catalog(), '2026_04')
        self.assertEqual(request['source'], BASE + 'Reference_Proteomes_2026_04.tar.gz')
        self.assertEqual((request['bytes'], request['md5']), (1234, '1' * 32))
        with self.assertRaises(ValueError):
            catalog_request(catalog(), '2026_03')
        with self.assertRaises(ValueError):
            catalog_request(catalog(digest='invalid'), '2026_04')

    def test_custom_source_is_verified_and_reused_without_an_archive_copy(self):
        source = self.root / 'input.faa'
        source.write_text('>protein OX=7955\nACDEFGHIKLMNPQRSTVWY\n')
        sha = hashlib.sha256(source.read_bytes()).hexdigest()
        request = resolve(self.root / 'db', 'fixture', str(source), sha, 'fasta')
        output = self.root / 'download'
        fetch(request, output)
        self.assertTrue((output / 'archive').is_symlink())
        record = validate_download(request, output / 'archive', output / 'source.json')
        self.assertEqual(record['sha256'], sha)
        with self.assertRaises(ValueError):
            validate_download(dict(request, sha256='0' * 64), output / 'archive', output / 'source.json')
        source.write_text('changed')
        with self.assertRaises(ValueError):
            validate_download(request, output / 'archive', output / 'source.json')

    def test_bad_source_checksum_never_publishes_success_metadata(self):
        source = self.root / 'input.faa'
        source.write_text('not the expected source')
        request = resolve(self.root / 'db', 'fixture', str(source), '0' * 64, 'fasta')
        output = self.root / 'bad'
        with self.assertRaises(ValueError):
            fetch(request, output)
        self.assertFalse((output / 'source.json').exists())

    def test_existing_database_reuses_offline_but_rejects_a_different_request(self):
        bundle = self.root / 'diamond/uniprot_reference_proteomes/2026_04'
        bundle.mkdir(parents=True)
        (bundle / 'reference_proteomes.dmnd').write_text('fixture index')
        identity = dict(collection='uniprot_reference_proteomes', release='2026_04',
                        source=BASE + 'Reference_Proteomes_2026_04.tar.gz',
                        source_format='tar', source_sha256='a' * 64)
        (bundle / 'protein_database.json').write_text(json.dumps(identity))
        with patch('contamination_protein_source.subprocess.check_output', side_effect=AssertionError('network')):
            self.assertTrue(resolve(self.root, '2026_04')['existing'])
        with self.assertRaises(ValueError):
            resolve(self.root, '2026_04', source_format='fasta')
        with self.assertRaises(ValueError):
            resolve(self.root, '2026_04', identity['source'], 'b' * 64)

    def test_catalog_resolution_and_custom_input_validation(self):
        request = resolve(self.root, '2026_04', catalog_text=catalog())
        self.assertFalse(request['existing'])
        self.assertEqual(request['md5'], '1' * 32)
        for release in ('../escape', '/absolute', ''):
            with self.assertRaises(ValueError):
                resolve(self.root, release)
        with self.assertRaises(ValueError):
            resolve(self.root, '2026_04', source='https://example.org/archive')

    def test_https_source_checks_the_publisher_md5_before_publishing(self):
        data = b'fixture archive bytes'
        request = catalog_request(catalog(digest=hashlib.md5(data).hexdigest()), '2026_04')
        request['bytes'] = len(data)
        def download(command):
            Path(command[command.index('--output') + 1]).write_bytes(data)
            return type('Completed', (), dict(returncode=0))()
        with patch('contamination_protein_source.subprocess.run', side_effect=download):
            output = self.root / 'verified'
            fetch(request, output)
            self.assertEqual((output / 'archive').read_bytes(), data)
            self.assertEqual(json.loads((output / 'source.json').read_text())['sha256'], hashlib.sha256(data).hexdigest())
            with self.assertRaises(ValueError):
                fetch(dict(request, md5='0' * 32), self.root / 'rejected')
            self.assertFalse((self.root / 'rejected/source.json').exists())

    def test_canonical_archive_filter_and_merged_taxids_preserve_sequences(self):
        sequence = 'ACDEFGHIKLMNPQRSTVWY' * 5
        archive = self.root / 'proteomes.tar.gz'
        with tarfile.open(archive, 'w:gz') as bundle:
            for name, text in [('species/UP000000001_10.fasta.gz', f'>canonical OX=10\n{sequence}\n'),
                               ('species/UP000000001_10_additional.fasta.gz', '>isoform without taxid\nAAAA\n'),
                               ('species/UP000000001_10_DNA.fasta.gz', '>CDS without taxid\nACGT\n')]:
                data = gzip.compress(text.encode())
                member = tarfile.TarInfo(name); member.size = len(data)
                bundle.addfile(member, io.BytesIO(data))
        nodes, merged = self.root / 'nodes.dmp', self.root / 'merged.dmp'
        nodes.write_text('1\t|\t1\t|\tno rank\t|\n20\t|\t1\t|\tspecies\t|\n')
        merged.write_text('10\t|\t20\t|\n')
        fasta, mapping, probe = [self.root / name for name in ('out.faa.gz', 'map.tsv.gz', 'probe.faa')]
        result = prepare(archive, 'tar', nodes, merged, fasta, mapping, probe)
        self.assertEqual(result, (1, 100, 1, 20))
        with gzip.open(fasta, 'rt') as handle:
            self.assertEqual(handle.read(), f'>UPR1 canonical OX=10\n{sequence}\n')
        with gzip.open(mapping, 'rt') as handle:
            self.assertEqual(handle.read(), 'accession\taccession.version\ttaxid\tgi\nUPR1\tUPR1\t20\t0\n')


if __name__ == '__main__':
    unittest.main()
