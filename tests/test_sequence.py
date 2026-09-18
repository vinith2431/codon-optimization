import unittest

from core.sequence import (
    clean_dna,
    validate_dna,
    get_codons,
    gc_content,
    reverse_complement,
    dna_to_rna
)


class TestSequence(unittest.TestCase):

    def test_clean_dna(self):
        sequence = " atg gcc "
        self.assertEqual(
            clean_dna(sequence),
            "ATGGCC"
        )

    def test_validate_dna(self):
        self.assertTrue(
            validate_dna("ATGGCC")
        )

    def test_invalid_dna(self):
        with self.assertRaises(ValueError):
            validate_dna("ATGXCC")

    def test_invalid_length(self):
        with self.assertRaises(ValueError):
            validate_dna("ATGG")

    def test_get_codons(self):
        self.assertEqual(
            get_codons("ATGGCCGCG"),
            ["ATG", "GCC", "GCG"]
        )

    def test_gc_content(self):
        self.assertAlmostEqual(
            gc_content("ATGGCC"),
            4 / 6
        )

    def test_reverse_complement(self):
        self.assertEqual(
            reverse_complement("ATGC"),
            "GCAT"
        )

    def test_dna_to_rna(self):
        self.assertEqual(
            dna_to_rna("ATGC"),
            "AUGC"
        )


if __name__ == "__main__":
    unittest.main()