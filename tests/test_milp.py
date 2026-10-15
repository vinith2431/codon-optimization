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

    def test_milp_forbidden_sites(self):
        # Ala-Ile-Phe-Glu: Ala can be GCA/GCC/GCG/GCT, Ile can be ATT/ATC/ATA, Phe can be TTT/TTC
        # If Ala=GAA... Glu is GAA/GAG.
        seq = "ATGGAATTCTTC"  # Met-Glu-Phe-Phe containing EcoRI "GAATTC"
        # Optimize forbidding GAATTC
        result = optimize_milp(
            seq,
            self.organism,
            forbidden_sites=["GAATTC"],
        )
        self.assertNotIn("GAATTC", result)
        self.assertTrue(protein_preserved(seq, result))

    def test_milp_homopolymer_avoidance(self):
        # Lys-Lys: AAA AAA -> 6 A's by default
        seq = "ATGAAAAAAAAA"  # Met-Lys-Lys-Lys
        result = optimize_milp(
            seq,
            self.organism,
            max_homopolymer=4,
        )
        self.assertNotIn("AAAAA", result)
        self.assertTrue(protein_preserved(seq, result))

    def test_milp_lazy_repeat_elimination(self):
        # 4 identical codons: Ala-Ala-Ala-Ala-Ala-Ala (18 bp)
        # Without repeat constraint, MILP selects the most preferred codon for all (repeating 18 bp)
        seq = "ATG" + "GCG" * 8  # 27 bp
        unconstrained = optimize_milp(seq, self.organism)
        # unconstrained will have repeated 12-mers
        from constraints.repeats import find_repeats
        self.assertTrue(len(find_repeats(unconstrained, 12)) > 0)

        # With repeat_length=12
        constrained = optimize_milp(seq, self.organism, repeat_length=12)
        repeats = find_repeats(constrained, 12)
        self.assertEqual(len(repeats), 0)
        self.assertTrue(protein_preserved(seq, constrained))


if __name__ == "__main__":
    unittest.main()