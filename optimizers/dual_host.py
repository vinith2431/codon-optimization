"""
Dual-host codon optimization.

The objective combines codon usage preferences from
two target organisms.
"""

from core.sequence import get_codons
from core.verification import CODON_TO_AA
from data.codon_tables import get_codon_table


def dual_host_score(
    chromosome,
    organism1,
    organism2,
    weight1=0.5,
    weight2=0.5
):
    """
    Calculate combined codon usage score for two hosts.
    """

    table1 = get_codon_table(
        organism1
    )

    table2 = get_codon_table(
        organism2
    )

    scores = []

    for codon in chromosome:

        amino_acid = CODON_TO_AA[codon]

        score1 = table1[
            amino_acid
        ].get(codon, 0.0)

        score2 = table2[
            amino_acid
        ].get(codon, 0.0)

        combined = (
            weight1 * score1
            + weight2 * score2
        )

        scores.append(combined)

    if not scores:
        return 0.0

    return sum(scores) / len(scores)


def optimize_dual_host(
    sequence,
    organism1,
    organism2,
    weight1=0.5,
    weight2=0.5
):
    """
    Perform a simple dual-host codon optimization.

    For every amino-acid position, select the synonymous
    codon with the highest weighted preference across
    the two hosts.
    """

    codons = get_codons(sequence)

    table1 = get_codon_table(
        organism1
    )

    table2 = get_codon_table(
        organism2
    )

    optimized = []

    for original_codon in codons:

        amino_acid = CODON_TO_AA[
            original_codon
        ]

        candidates = set(
            table1[amino_acid].keys()
        ) | set(
            table2[amino_acid].keys()
        )

        best_codon = None
        best_score = float("-inf")

        for codon in candidates:

            score1 = table1[
                amino_acid
            ].get(codon, 0.0)

            score2 = table2[
                amino_acid
            ].get(codon, 0.0)

            score = (
                weight1 * score1
                + weight2 * score2
            )

            if score > best_score:

                best_score = score
                best_codon = codon

        optimized.append(
            best_codon
        )

    return "".join(optimized)