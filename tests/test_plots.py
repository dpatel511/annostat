from __future__ import annotations

import unittest

import matplotlib as mpl
from matplotlib.backends.backend_agg import FigureCanvasAgg

from annostat.plots import _STYLE, _new_figure


class PlotLayoutTests(unittest.TestCase):
    def test_header_spacing_is_safe_for_short_and_tall_figures(self) -> None:
        for height in (3.2, 4.32, 6.5, 11.5):
            for dpi in (72, 100, 200):
                with self.subTest(height=height, dpi=dpi), mpl.rc_context(_STYLE):
                    figure = _new_figure(
                        12, height, "Start codon usage",
                        "Observed first codon across 4,237 coding sequences",
                    )
                    figure.set_dpi(dpi)
                    axis = figure.subplots()
                    canvas = FigureCanvasAgg(figure)
                    canvas.draw()
                    renderer = canvas.get_renderer()
                    title = figure._suptitle.get_window_extent(renderer)
                    description = figure.texts[1].get_window_extent(renderer)
                    axes = axis.get_window_extent(renderer)
                    self.assertGreaterEqual(title.y0 - description.y1, 3 * dpi / 72)
                    self.assertGreaterEqual(description.y0 - axes.y1, 3 * dpi / 72)
                    self.assertLessEqual(title.y1, figure.bbox.y1)
                    self.assertGreaterEqual(figure.bbox.y1 - title.y1, 3 * dpi / 72)
                    self.assertGreaterEqual(description.y0, figure.bbox.y0)
                    figure.clear()


if __name__ == "__main__":
    unittest.main()
