"""Tests for exact fraction parsing and label formatting."""

import importlib.util
import unittest
from decimal import Decimal
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "bin" / "busco_multigene_tree.py"
SPEC = importlib.util.spec_from_file_location("busco_multigene_tree", MODULE_PATH)
busco_multigene_tree = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(busco_multigene_tree)


class FractionLabelTests(unittest.TestCase):
    """Exercise exact decimal fraction label behaviour."""

    def test_decimal_fraction_uses_underscore_label(self) -> None:
        """Decimal percentages use underscores in output labels."""
        self.assertEqual(
            busco_multigene_tree.fraction_to_label(Decimal("0.999")),
            "frac99_9pct",
        )

    def test_whole_fraction_omits_decimal_separator(self) -> None:
        """Whole percentages keep the compact integer label."""
        self.assertEqual(
            busco_multigene_tree.fraction_to_label(Decimal("1.0")),
            "frac100pct",
        )

    def test_parse_fractions_deduplicates_and_sorts(self) -> None:
        """Parsed fractions are normalised, deduplicated, and sorted."""
        self.assertEqual(
            busco_multigene_tree.parse_fractions("1.0,0.999,0.80,0.8"),
            [Decimal("0.8"), Decimal("0.999"), Decimal("1")],
        )


if __name__ == "__main__":
    unittest.main()
