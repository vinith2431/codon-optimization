"""
MILP-based codon optimization.

The optimizer selects one synonymous codon for every
amino-acid position while maximizing codon usage preference.

Global GC-content constraints can be applied.
"""

import pulp

from core.sequence import get_codons
from core.verification import CODON_TO_AA
from data.codon_tables import get_codon_table


def optimize_milp(
    sequence,
    organism,
    gc_min=None,
    gc_max=None
):
    """
    Optimize a DNA sequence using MILP.

    Parameters
    ----------
    sequence : str
        Original DNA sequence.

    organism : str
        Target organism.

    gc_min : float, optional
        Minimum global GC fraction.

    gc_max : float, optional
        Maximum global GC fraction.

    Returns
    -------
    str
        Optimized DNA sequence.
    """

    original_codons = get_codons(sequence)

    table = get_codon_table(organism)

    problem = pulp.LpProblem(
        "Codon_Optimization",
        pulp.LpMaximize
    )

    variables = {}

    # --------------------------------------------------
    # Create binary decision variables
    # --------------------------------------------------

    for position, original_codon in enumerate(original_codons):

        amino_acid = CODON_TO_AA[original_codon]

        synonymous_codons = table[amino_acid]

        variables[position] = {}

        for codon in synonymous_codons:

            variables[position][codon] = pulp.LpVariable(
                f"x_{position}_{codon}",
                cat=pulp.LpBinary
            )

    # --------------------------------------------------
    # Constraint:
    # Exactly one codon must be selected per position
    # --------------------------------------------------

    for position in variables:

        problem += (
            pulp.lpSum(
                variables[position].values()
            ) == 1
        )

    # --------------------------------------------------
    # Objective:
    # Maximize codon usage preference
    # --------------------------------------------------

    objective_terms = []

    for position in variables:

        for codon, variable in variables[position].items():

            frequency = table[
                CODON_TO_AA[
                    codon
                ]
            ][codon]

            objective_terms.append(
                frequency * variable
            )

    problem += pulp.lpSum(objective_terms)

    # --------------------------------------------------
    # Global GC-content constraint
    # --------------------------------------------------

    total_bases = len(sequence)

    gc_expression = []

    for position in variables:

        for codon, variable in variables[position].items():

            gc_count = (
                codon.count("G")
                + codon.count("C")
            )

            gc_expression.append(
                gc_count * variable
            )

    total_gc = pulp.lpSum(gc_expression)

    if gc_min is not None:

        problem += (
            total_gc >= gc_min * total_bases
        )

    if gc_max is not None:

        problem += (
            total_gc <= gc_max * total_bases
        )

    # --------------------------------------------------
    # Solve
    # --------------------------------------------------

    solver = pulp.PULP_CBC_CMD(
        msg=False
    )

    status = problem.solve(solver)

    if pulp.LpStatus[status] != "Optimal":

        raise RuntimeError(
            "MILP optimization did not find an optimal solution. "
            f"Solver status: {pulp.LpStatus[status]}"
        )

    # --------------------------------------------------
    # Recover optimized sequence
    # --------------------------------------------------

    optimized_codons = []

    for position in range(len(original_codons)):

        selected_codon = None

        for codon, variable in variables[position].items():

            if pulp.value(variable) > 0.5:

                selected_codon = codon
                break

        if selected_codon is None:

            raise RuntimeError(
                f"No codon selected at position {position}."
            )

        optimized_codons.append(
            selected_codon
        )

    return "".join(optimized_codons)


def milp_optimization(
    sequence,
    organism,
    gc_min=None,
    gc_max=None
):
    """
    Alias for optimize_milp().
    """

    return optimize_milp(
        sequence=sequence,
        organism=organism,
        gc_min=gc_min,
        gc_max=gc_max
    )