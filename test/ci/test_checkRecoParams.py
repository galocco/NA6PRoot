"""Exercise static validation without ROOT or a project build."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from checkRecoParams import REPOSITORY, read_schema, validate_file


class ParameterValidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schema = read_schema(REPOSITORY / 'rec/include/NA6PRecoParam.h')

    def check_config(self, contents):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'params.ini'
            path.write_text(contents)
            return validate_file(path, self.schema)

    def test_valid_partial_configuration(self):
        self.assertEqual(self.check_config(
            '[reco]\nvtNLayers=5\nuseLinRefVT=true\n'
            'maxPropagationStep=2e-3\nvertexerRecoType=YZ\n'
            'vtMaxDeltaPhiCellsCA[9]=999\nmtNSigma2Crude[4]=9\n'), [])

    def test_invalid_values_and_names(self):
        for assignment in (
                'msZForMSMIDmatch=795', 'VTNLayers=5', 'vtNLayers=2.5',
                'vtNLayers=2147483648', 'useLinRefVT=maybe',
                'maxPropagationStep=nan', 'maxPropagationStep=1e100',
                'vtNLayers[0]=5', 'mtNSigma2Crude=9',
                'mtNSigma2Crude[5]=9', 'vtMinNClusTracksCA[10]=3',
                'vtMinNClusTracksCA[-1]=3', 'vertexerRecoType='):
            with self.subTest(assignment=assignment):
                self.assertTrue(self.check_config('[reco]\n' + assignment + '\n'))

    def test_malformed_or_empty_configuration(self):
        for contents in ('', '[reco]\n', 'vtNLayers=5\n',
                         '[other]\nvtNLayers=5\n',
                         '[reco]\nvtNLayers=5\nvtNLayers=4\n',
                         '[DEFAULT]\nvtNLayers=5\n[reco]\n'):
            with self.subTest(contents=contents):
                self.assertTrue(self.check_config(contents))

    def test_cli_discovers_nested_files_and_returns_failure(self):
        script = Path(__file__).with_name('checkRecoParams.py')
        with tempfile.TemporaryDirectory() as directory:
            command = [sys.executable, str(script), '--params', directory]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('No .ini files found', result.stderr)
            nested = Path(directory) / 'nested'
            nested.mkdir()
            config = nested / 'example.ini'
            config.write_text('[reco]\nvtNLayers=5\n')
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('1/1 parameter files passed', result.stdout)
            config.write_text('[reco]\nunknownParameter=1\n')
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('unknownParameter: unknown parameter', result.stderr)


if __name__ == '__main__':
    unittest.main()
