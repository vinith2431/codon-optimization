import unittest

from constraints.gc import (
    calculate_gc,
    gc_in_range,
    gc_constraint_violation
)

from constraints.local_gc import (
    local_gc_content,
    local_gc_in_range
)

from constraints.repeats import (
    find_repeats,
    has_repeats
)

from constraints.homopolymer import (
    find_homopolymers,
    has_homopolymer
)

from constraints.restriction_sites import (
    find_restriction_sites,
    contains_restriction_site
)


class TestConstraints(unittest.TestCase):

    # ----------------------------------------------
    # Global GC
    # ----------------------------------------------

    def test_gc(self):

        self.assertAlmostEqual(
            calculate_gc("ATGGCC"),
            4 / 6
        )

    def test_gc_range_valid(self):

        self.assertTrue(
            gc_in_range(
                "ATGGCC",
                0.5,
                0.8
            )
        )

    def test_gc_range_invalid(self):

        self.assertFalse(
            gc_in_range(
                "AAAAAA",
                0.5,
                0.8
            )
        )

    def test_gc_violation(self):

        violation = gc_constraint_violation(
            "AAAAAA",
            0.5,
            0.8
        )

        self.assertAlmostEqual(
            violation,
            0.5
        )

    # ----------------------------------------------
    # Local GC
    # ----------------------------------------------

    def test_local_gc(self):

        result = local_gc_content(
            "ATGGCC",
            3
        )

        self.assertEqual(
            len(result),
            4
        )

    def test_local_gc_range(self):

        self.assertTrue(
            local_gc_in_range(
                "ATGGCC",
                3,
                0.0,
                1.0
            )
        )

    # ----------------------------------------------
    # Repeats
    # ----------------------------------------------

    def test_repeats(self):

        result = find_repeats(
            "ATGCATGC",
            4
        )

        self.assertTrue(
            len(result) > 0
        )

    def test_has_repeats(self):

        self.assertTrue(
            has_repeats(
                "ATGCATGC",
                4
            )
        )

    # ----------------------------------------------
    # Homopolymers
    # ----------------------------------------------

    def test_homopolymer(self):

        result = find_homopolymers(
            "ATGAAAAAA",
            5
        )

        self.assertTrue(
            len(result) > 0
        )

    def test_has_homopolymer(self):

        self.assertTrue(
            has_homopolymer(
                "ATGAAAAAA",
                5
            )
        )

    # ----------------------------------------------
    # Restriction sites
    # ----------------------------------------------

    def test_restriction_site(self):

        result = find_restriction_sites(
            "ATGGAATTCGCC",
            ["GAATTC"]
        )

        self.assertEqual(
            len(result),
            1
        )

    def test_contains_restriction_site(self):

        self.assertTrue(
            contains_restriction_site(
                "ATGGAATTCGCC",
                ["GAATTC"]
            )
        )


if __name__ == "__main__":
    unittest.main()