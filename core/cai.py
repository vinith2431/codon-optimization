"""
Codon Adaptation Index (CAI) calculation.
"""

import math

from core.sequence import get_codons
from data.codon_tables import get_codon_table


def calculate_cai(sequence, organism):
    """
    Calculate the Codon Adaptation Index (CAI).

    CAI is calculated using the geometric mean of
    relative synonymous codon usage values.
    """

    codons = get_codons(sequence)

    table = get_codon_table(organism)

    relative_values = []

    for codon in codons:

        codon = codon.upper()

        amino_acid = None

        for aa, synonymous in table.items():

            if codon in synonymous:
                amino_acid = aa
                break

        if amino_acid is None:
            raise ValueError(
                f"Codon '{codon}' not found."
            )

        synonymous_codons = table[amino_acid]

        maximum_frequency = max(
            synonymous_codons.values()
        )

        codon_frequency = synonymous_codons[codon]

        relative_value = (
            codon_frequency / maximum_frequency
        )

        relative_values.append(relative_value)

    if not relative_values:
        return 0.0

    # Avoid log(0)
    epsilon = 1e-10

    log_sum = sum(
        math.log(max(value, epsilon))
        for value in relative_values
    )

    cai = math.exp(
        log_sum / len(relative_values)
    )

    return cai