import unittest

from optimizers.greedy import (
    optimize_greedy
)

from core.verification import (
    protein_preserved
)


class TestGreedy(unittest.TestCase):

    def setUp(self):

        self.sequence = (
            "ATGGCCGCG"
        )

        self.organism = "ecoli"

    def test_optimization_returns_sequence(self):

        result = optimize_greedy(
            self.sequence,
            self.organism
        )

        self.assertIsInstance(
            result,
            str
        )

    def test_same_length(self):

        result = optimize_greedy(
            self.sequence,
            self.organism
        )

        self.assertEqual(
            len(result),
            len(self.sequence)
        )

    def test_protein_preserved(self):

        result = optimize_greedy(
            self.sequence,
            self.organism
        )

        self.assertTrue(
            protein_preserved(
                self.sequence,
                result
            )
        )

    def test_valid_dna(self):

        result = optimize_greedy(
            self.sequence,
            self.organism
        )

        self.assertTrue(
            all(
                base in "ATGC"
                for base in result
            )
        )


if __name__ == "__main__":
    unittest.main()