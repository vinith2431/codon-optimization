import unittest

from optimizers.milp import (
    optimize_milp
)

from core.verification import (
    protein_preserved
)

from constraints.gc import (
    gc_in_range
)


class TestMILP(unittest.TestCase):

    def setUp(self):

        self.sequence = (
            "ATGGCCGCG"
        )

        self.organism = "ecoli"

    def test_milp_returns_sequence(self):

        result = optimize_milp(
            self.sequence,
            self.organism
        )

        self.assertIsInstance(
            result,
            str
        )

    def test_same_length(self):

        result = optimize_milp(
            self.sequence,
            self.organism
        )

        self.assertEqual(
            len(result),
            len(self.sequence)
        )

    def test_protein_preserved(self):

        result = optimize_milp(
            self.sequence,
            self.organism
        )

        self.assertTrue(
            protein_preserved(
                self.sequence,
                result
            )
        )

    def test_gc_constraint(self):

        result = optimize_milp(
            self.sequence,
            self.organism,
            gc_min=0.40,
            gc_max=0.70
        )

        self.assertTrue(
            gc_in_range(
                result,
                0.40,
                0.70
            )
        )


if __name__ == "__main__":
    unittest.main()