"""Verify CLI stream separation, mode disclosure, and expected failures."""
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from annostat.cli import main


class CliTests(unittest.TestCase):
    def invoke(self, arguments):
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            try:
                code = main(arguments)
            except SystemExit as error:
                code = error.code
        return code, stdout.getvalue(), stderr.getvalue()

    def compare_arguments(self, *extra):
        return ['compare', '--genome', 'A', 'genome.fna', 'a.gff3',
                '--genome', 'B', 'genome.fna', 'b.gff3', *extra]

    def comparison(self, gene=False):
        result = {'datasets': [{}, {}], 'pairwise_comparisons': [{}],
                  'taxonomic_scope': 'not_determined', 'warnings': ['Source context differs.']}
        if gene:
            result['gene_level_comparison'] = {
                'relationships': [{}, {}, {}], 'exact_model_count': 1,
                'status_counts': {'matched': 2, 'unmatched': 1},
                'consequence_counts': {'identical': 1, 'terminal_only': 1},
                'cds_counts': {'A': 3, 'B': 2}, 'unmatched_cds_counts': {'A': 1, 'B': 0},
            }
        return result

    def test_aggregate_mode_is_explicit_and_progress_uses_stderr(self):
        def compare(datasets, output, progress, **kwargs):
            self.assertFalse(kwargs['gene_level'])
            progress('Checking inputs')
            return self.comparison()

        with patch('annostat.cli.run_comparison', side_effect=compare):
            code, stdout, stderr = self.invoke(self.compare_arguments())
        self.assertEqual(code, 0)
        self.assertIn('aggregate profiles (no CDS-by-CDS matching)', stderr)
        self.assertIn('[1/4] Checking inputs', stderr)
        self.assertNotIn('[1/4]', stdout)
        self.assertIn('Individual CDS models were not matched', stdout)
        self.assertIn('note: Source context differs.', stderr)

    def test_gene_summary_exposes_denominators_and_per_source_counts(self):
        with patch('annostat.cli.run_comparison', return_value=self.comparison(gene=True)) as compare:
            code, stdout, stderr = self.invoke(self.compare_arguments('--gene-level'))
        self.assertEqual(code, 0)
        self.assertTrue(compare.call_args.kwargs['gene_level'])
        self.assertIn('gene-level (CDS correspondence and protein evidence)', stderr)
        self.assertIn('Accepted CDS pairs', stdout)
        self.assertIn('1 (50.00% of accepted pairs)', stdout)
        self.assertIn('A: 3 CDS; 1 without an accepted pair; 1 unmatched', stdout)
        self.assertIn('comparison.txt', stdout)
        self.assertIn('Agreement is not annotation accuracy', stdout)

    def test_empty_gene_summary_does_not_report_zero_as_measured_percentage(self):
        summary = self.comparison(gene=True)
        gene = summary['gene_level_comparison']
        gene.update(relationships=[], exact_model_count=0, status_counts={}, consequence_counts={},
                    cds_counts={'A': 0, 'B': 0}, unmatched_cds_counts={'A': 0, 'B': 0})
        with patch('annostat.cli.run_comparison', return_value=summary):
            code, stdout, stderr = self.invoke(self.compare_arguments('--gene-level'))
        self.assertEqual(code, 0)
        self.assertIn('not available', stdout)
        self.assertNotIn('0.00%', stdout)

    def test_quiet_compare_suppresses_success_output_but_retains_errors(self):
        with patch('annostat.cli.run_comparison', return_value=self.comparison()):
            code, stdout, stderr = self.invoke(self.compare_arguments('--quiet'))
        self.assertEqual((code, stdout, stderr), (0, '', ''))
        with patch('annostat.cli.run_comparison', side_effect=ValueError('genomes differ')):
            code, stdout, stderr = self.invoke(self.compare_arguments('--quiet'))
        self.assertEqual(code, 2)
        self.assertEqual(stdout, '')
        self.assertIn('annostat compare: error: genomes differ', stderr)
        self.assertNotIn('usage:', stderr)
        self.assertNotIn('Traceback', stderr)

    def test_argument_errors_keep_usage_and_explain_pairwise_requirement(self):
        with patch('annostat.cli.run_comparison') as compare:
            code, stdout, stderr = self.invoke(self.compare_arguments(
                '--genome', 'C', 'genome.fna', 'c.gff3', '--gene-level'))
        self.assertEqual(code, 2)
        self.assertIn('usage:', stderr)
        self.assertIn('run each pair separately', stderr)
        compare.assert_not_called()

    def test_runtime_errors_are_consistent_across_commands(self):
        for command, target in (
            (['analyze', '-f', 'a.fna', '-g', 'a.gff3', '-q'], 'run_analysis'),
            (['summarize', 'results', '-q'], 'build_cohort'),
            (['fetch', 'GCF_000009045.1'], 'fetch_genomes'),
        ):
            with self.subTest(command=command[0]), patch(
                f'annostat.cli.{target}', side_effect=PermissionError(13, 'Permission denied', 'blocked')
            ):
                code, stdout, stderr = self.invoke(command)
                self.assertEqual(code, 2)
                self.assertEqual(stdout, '')
                self.assertIn("cannot access 'blocked': Permission denied", stderr)
                self.assertNotIn('usage:', stderr)
                self.assertNotIn('Traceback', stderr)

    def test_validation_output_failure_is_a_concise_runtime_error(self):
        with patch('annostat.cli.validate_annotation', return_value={
            'severity_counts': {'error': 1, 'warning': 0, 'info': 0},
        }), patch('annostat.cli.write_validation', side_effect=OSError('cannot write results')):
            code, stdout, stderr = self.invoke(['validate', '-f', 'a', '-g', 'b'])
        self.assertEqual(code, 2)
        self.assertIn('annostat validate: error: cannot write results', stderr)
        self.assertNotIn('usage:', stderr)

    def test_validation_warning_threshold_and_quiet_error_diagnostics(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta, gff = root / 'a.fna', root / 'a.gff3'
            fasta.write_text('>chr\nATGAAATAA\n', encoding='utf-8')
            gff.write_text('chr\ttest\tCDS\t1\t9\t.\t+\t0\tID=c\n', encoding='utf-8')
            arguments = ['validate', '-f', str(fasta), '-g', str(gff), '-o', str(root / 'out')]
            code, stdout, stderr = self.invoke(arguments + ['--fail-on', 'warning'])
            self.assertEqual(code, 1)
            self.assertIn('threshold reached', stdout)
            self.assertIn('GFF_VERSION', stderr)
            code, stdout, stderr = self.invoke(arguments + ['--quiet'])
            self.assertEqual((code, stdout, stderr), (0, '', ''))
            gff.write_text('bad\n', encoding='utf-8')
            code, stdout, stderr = self.invoke(arguments + ['--quiet'])
            self.assertEqual(code, 1)
            self.assertEqual(stdout, '')
            self.assertIn('GFF_PARSE', stderr)
            self.assertIn('threshold reached', stderr)

    def test_validation_preview_is_bounded_and_prioritizes_errors(self):
        findings = [{'severity': 'warning', 'rule_id': 'WARN', 'message': 'review context'}]
        findings += [{'severity': 'error', 'rule_id': f'ERR{i}', 'message': 'invalid record'}
                     for i in range(7)]
        result = {'valid': False, 'severity_counts': {'error': 7, 'warning': 1, 'info': 0},
                  'findings': findings}
        with patch('annostat.cli.validate_annotation', return_value=result), patch(
            'annostat.cli.write_validation', return_value=[Path('validation.json'), Path('validation.tsv')]
        ):
            code, stdout, stderr = self.invoke(['validate', '-f', 'a', '-g', 'b'])
        self.assertEqual(code, 1)
        self.assertEqual(stderr.count('invalid record'), 5)
        self.assertNotIn('review context', stderr)
        self.assertIn('3 more diagnostic(s)', stderr)

    def test_fetch_quiet_suppresses_success_summary(self):
        with patch('annostat.cli.fetch_genomes', return_value=[]):
            code, stdout, stderr = self.invoke(['fetch', 'GCF_000009045.1', '--quiet'])
        self.assertEqual((code, stdout, stderr), (0, '', ''))


if __name__ == '__main__':
    unittest.main()
