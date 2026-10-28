"""Check plotted measurements against small independent, known-value examples."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from annostat import plots


class PlotDataTests(unittest.TestCase):
    def render(self, function, check, *args, **kwargs):
        original = plots._save_svg

        def inspect(figure, path, title, description):
            check(figure)
            original(figure, path, title, description)

        with tempfile.TemporaryDirectory() as temporary, patch('annostat.plots._save_svg', side_effect=inspect):
            path = Path(temporary) / 'plot.svg'
            function(path, *args, **kwargs)
            self.assertGreater(path.stat().st_size, 0)

    def test_gene_bars_use_accepted_pairs_and_preserve_zero_categories(self):
        result = {'status_counts': {'matched': 4, 'unmatched': 3},
                  'consequence_counts': {'identical': 2, 'terminal_only': 1, 'limited': 1}}

        def check(figure):
            axis = figure.axes[0]
            self.assertEqual([bar.get_width() for bar in axis.patches], [2, 1, 0, 0, 1])
            self.assertEqual([text.get_text() for text in axis.texts], ['2 / 4', '1 / 4', '0 / 4', '0 / 4', '1 / 4'])

        self.render(plots.write_gene_consequences, check, result)

    def test_bar_labels_use_the_supplied_cds_denominator(self):
        def check(figure):
            axis = figure.axes[0]
            self.assertEqual([bar.get_width() for bar in axis.patches], [8, 5])
            self.assertEqual([text.get_text() for text in axis.texts], ['8 (80.00%)', '5 (50.00%)'])

        self.render(plots.write_bar_chart, check, 'COG categories', {'J': 5, 'K': 8}, percentage_total=10)

    def test_histogram_cutoff_does_not_change_full_sample_summary(self):
        values = [3] * 99 + [3000]

        def check(figure):
            axis = figure.axes[0]
            self.assertEqual(sum(bar.get_height() for bar in axis.patches), 99)
            self.assertEqual(list(axis.lines[0].get_xdata()), [3, 3])
            self.assertIn('Mean: 33 nt (outside view)', [line.get_label() for line in axis.lines])
            self.assertIn('1 CDS above 3 nt', figure.texts[1].get_text())

        self.render(plots.write_histogram, check, 'CDS lengths', values)

    def test_overview_keeps_measured_zero_and_missing_values_separate(self):
        keys = ('gc_percent', 'coding_density_percent', 'cog_coverage_percent',
                'hypothetical_percent', 'gene_name_percent', 'cds_per_mb')
        profiles = [{'label': label, **dict.fromkeys(keys, value)}
                    for label, value in (('zero', 0), ('missing', None), ('measured', 50))]

        def check(figure):
            self.assertEqual(len(figure.axes), 6)
            for axis in figure.axes:
                self.assertEqual([list(point.get_offsets()[0]) for point in axis.collections], [[0, 0], [50, 2]])
                self.assertIn('Not available', [text.get_text() for text in axis.texts])

        self.render(plots.write_comparison_overview, check, profiles)

    def test_cog_deltas_are_percentage_points_b_minus_a(self):
        profiles = [{'label': 'A', 'cog_category_percentages': {'J': 25, 'K': 75}},
                    {'label': 'B', 'cog_category_percentages': {'J': 100, 'K': 0}}]

        def check(figure):
            self.assertEqual([bar.get_width() for bar in figure.axes[0].patches], [75, -75])
            self.assertEqual([text.get_text() for text in figure.axes[0].texts], ['+75.00 pp', '-75.00 pp'])

        self.render(plots.write_cog_comparison, check, profiles, {'J': 'Translation', 'K': 'Transcription'})

    def test_heatmap_values_and_text_contrast_follow_the_actual_colors(self):
        profiles = [{'label': label, 'cog_category_percentages': {'J': value, 'K': 100 - value}}
                    for label, value in (('A', 0), ('B', 25), ('C', 100))]

        def check(figure):
            axis = figure.axes[0]
            np.testing.assert_array_equal(axis.images[0].get_array(), [[0, 25, 100], [100, 75, 0]])
            for text in axis.texts:
                if text.get_text() == '0.00%':
                    self.assertEqual(text.get_color(), 'white')
                if text.get_text() == '100.00%':
                    self.assertEqual(text.get_color(), plots._INK)

        self.render(plots.write_cog_comparison, check, profiles, {'J': 'Translation', 'K': 'Transcription'})


if __name__ == '__main__':
    unittest.main()
