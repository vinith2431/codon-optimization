"""
Greedy codon optimization.

For each amino acid, the codon with the highest usage
frequency is selected.
"""

from core.sequence import get_codons
from core.verification import CODON_TO_AA
from data.codon_tables import get_codon_table


def optimize_greedy(sequence, organism):
    """
    Optimize a DNA sequence using the highest-frequency
    synonymous codon for each amino acid.

    Parameters
    ----------
    sequence : str
        Original DNA sequence.

    organism : str
        Target organism.

    Returns
    -------
    str
        Greedy optimized DNA sequence.
    """

    original_codons = get_codons(sequence)

    table = get_codon_table(organism)

    optimized_codons = []

    for codon in original_codons:

        codon = codon.upper()

        if codon not in CODON_TO_AA:
            raise ValueError(
                f"Unknown codon: {codon}"
            )

        amino_acid = CODON_TO_AA[codon]

        # Preserve stop codons
        if amino_acid == "*":
            optimized_codons.append(codon)
            continue

        synonymous = table[amino_acid]

        best_codon = max(
            synonymous,
            key=synonymous.get
        )

        optimized_codons.append(best_codon)

    return "".join(optimized_codons)


def greedy_optimization(sequence, organism):
    """
    Alias for optimize_greedy().
    """

    return optimize_greedy(
        sequence,
        organism
    )