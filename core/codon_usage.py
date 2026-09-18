"""
Utilities for working with codon usage frequencies.
"""

from data.codon_tables import get_codon_table


def codon_frequency(codon, organism):
    """
    Return the usage frequency of a codon in an organism.
    """
    codon = codon.upper()

    table = get_codon_table(organism)

    for amino_acid, codons in table.items():

        if codon in codons:
            return codons[codon]

    raise ValueError(
        f"Codon '{codon}' not found in organism '{organism}'."
    )


def preferred_codon(amino_acid, organism):
    """
    Return the most frequently used codon
    for a given amino acid.
    """
    table = get_codon_table(organism)

    amino_acid = amino_acid.upper()

    if amino_acid not in table:
        raise ValueError(
            f"Unknown amino acid: {amino_acid}"
        )

    codons = table[amino_acid]

    return max(
        codons,
        key=codons.get
    )


def synonymous_codons(amino_acid, organism):
    """
    Return all synonymous codons for an amino acid.
    """
    table = get_codon_table(organism)

    amino_acid = amino_acid.upper()

    if amino_acid not in table:
        raise ValueError(
            f"Unknown amino acid: {amino_acid}"
        )

    return table[amino_acid]


def codon_usage_score(sequence, organism):
    """
    Calculate the average codon usage frequency
    of a DNA sequence.
    """

    from core.sequence import get_codons

    codons = get_codons(sequence)

    frequencies = [
        codon_frequency(codon, organism)
        for codon in codons
    ]

    if not frequencies:
        return 0.0

    return sum(frequencies) / len(frequencies)


def get_usage_table(organism):
    """
    Return the complete codon usage table.
    """
    return get_codon_table(organism)