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
from core.validator import assess_candidate


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

    def test_infeasible_constraints_raises_runtime_error(self):
        # ATGGCCGCG can never have GC >= 0.95 (Met is ATG, Ala can be GCT/GCC/GCA/GCG with 2 GC max each)
        with self.assertRaises(RuntimeError) as ctx:
            optimize_ga(
                self.sequence,
                self.organism,
                population_size=10,
                generations=5,
                gc_min=0.95,
                seed=42,
            )
        self.assertIn("satisfying all configured constraints", str(ctx.exception))

    def test_returned_candidate_is_feasible_by_assess_candidate(self):
        result = optimize_ga(
            self.sequence,
            self.organism,
            population_size=20,
            generations=10,
            gc_min=0.40,
            gc_max=0.80,
            max_homopolymer=4,
            seed=42,
        )
        res = assess_candidate(
            self.sequence,
            result,
            {
                "gc_min": 0.40,
                "gc_max": 0.80,
                "local_min": None,
                "local_max": None,
                "five_prime_min": None,
                "five_prime_max": None,
                "max_homopolymer": 4,
            },
        )
        self.assertTrue(res["feasible"])


if __name__ == "__main__":
    unittest.main()