import unittest

from optimizers.genetic_algorithm import (
    optimize_ga
)

from core.verification import (
    protein_preserved
)

from constraints.gc import (
    gc_in_range
)


class TestGA(unittest.TestCase):

    def setUp(self):

        self.sequence = (
            "ATGGCCGCG"
        )

        self.organism = "ecoli"

    def test_ga_returns_sequence(self):

        result = optimize_ga(
            self.sequence,
            self.organism,
            population_size=10,
            generations=5,
            seed=42
        )

        self.assertIsInstance(
            result,
            str
        )

    def test_same_length(self):

        result = optimize_ga(
            self.sequence,
            self.organism,
            population_size=10,
            generations=5,
            seed=42
        )

        self.assertEqual(
            len(result),
            len(self.sequence)
        )

    def test_protein_preserved(self):

        result = optimize_ga(
            self.sequence,
            self.organism,
            population_size=10,
            generations=5,
            seed=42
        )

        self.assertTrue(
            protein_preserved(
                self.sequence,
                result
            )
        )

    def test_gc_constraint(self):

        result = optimize_ga(
            self.sequence,
            self.organism,
            population_size=10,
            generations=5,
            gc_min=0.40,
            gc_max=0.70,
            seed=42
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