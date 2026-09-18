"""
Translation and optimization verification utilities.
"""

from core.sequence import (
    clean_dna,
    validate_dna,
    get_codons
)


CODON_TO_AA = {
    # Phenylalanine
    "TTT": "F",
    "TTC": "F",

    # Leucine
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",

    # Isoleucine
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",

    # Methionine
    "ATG": "M",

    # Valine
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",

    # Serine
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",

    # Proline
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",

    # Threonine
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",

    # Alanine
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",

    # Tyrosine
    "TAT": "Y",
    "TAC": "Y",

    # Histidine
    "CAT": "H",
    "CAC": "H",

    # Glutamine
    "CAA": "Q",
    "CAG": "Q",

    # Asparagine
    "AAT": "N",
    "AAC": "N",

    # Lysine
    "AAA": "K",
    "AAG": "K",

    # Aspartic acid
    "GAT": "D",
    "GAC": "D",

    # Glutamic acid
    "GAA": "E",
    "GAG": "E",

    # Cysteine
    "TGT": "C",
    "TGC": "C",

    # Tryptophan
    "TGG": "W",

    # Arginine
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",

    # Glycine
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",

    # Stop codons
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
}


def translate_dna(sequence, stop_at_stop=True):
    """
    Translate a DNA sequence into an amino-acid sequence.
    """

    sequence = clean_dna(sequence)

    validate_dna(sequence)

    codons = get_codons(sequence)

    protein = []

    for codon in codons:

        if codon not in CODON_TO_AA:
            raise ValueError(
                f"Unknown codon: {codon}"
            )

        amino_acid = CODON_TO_AA[codon]

        if amino_acid == "*" and stop_at_stop:
            break

        protein.append(amino_acid)

    return "".join(protein)


def verify_translation(original_sequence, optimized_sequence):
    """
    Verify that two DNA sequences encode the same protein.
    """

    original_protein = translate_dna(
        original_sequence
    )

    optimized_protein = translate_dna(
        optimized_sequence
    )

    return {
        "valid": original_protein == optimized_protein,
        "original_protein": original_protein,
        "optimized_protein": optimized_protein,
    }


def protein_preserved(original_sequence, optimized_sequence):
    """
    Return True if optimization preserves the protein sequence.
    """

    result = verify_translation(
        original_sequence,
        optimized_sequence
    )

    return result["valid"]