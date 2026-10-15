"""
MILP-based codon optimization.

The optimizer selects one synonymous codon for every
amino-acid position while maximizing codon usage preference.

Supported linear constraints:
- Global GC-content bounds
- 5' region GC-content bounds
- Local sliding-window GC-content bounds
- Forbidden restriction site elimination cuts
- Excessive homopolymer run elimination cuts
"""

import itertools
import pulp

from core.sequence import get_codons, clean_dna
from core.verification import CODON_TO_AA
from data.codon_tables import get_codon_table


def optimize_milp(
    sequence,
    organism,
    gc_min=None,
    gc_max=None,
    five_prime_gc_min=None,
    five_prime_gc_max=None,
    five_prime_window=60,
    local_gc_min=None,
    local_gc_max=None,
    local_gc_window=60,
    forbidden_sites=None,
    max_homopolymer=None,
):
    """
    Optimize a DNA sequence using Mixed Integer Linear Programming (MILP).

    Parameters
    ----------
    sequence : str
        Original coding DNA sequence.
    organism : str
        Target organism.
    gc_min : float, optional
        Minimum global GC fraction (e.g. 0.40).
    gc_max : float, optional
        Maximum global GC fraction (e.g. 0.60).
    five_prime_gc_min : float, optional
        Minimum GC fraction in the 5' region (e.g. 0.40).
    five_prime_gc_max : float, optional
        Maximum GC fraction in the 5' region (e.g. 0.50).
    five_prime_window : int, optional
        Length in bp of the 5' window (default 60).
    local_gc_min : float, optional
        Minimum GC fraction in local windows (e.g. 0.30).
    local_gc_max : float, optional
        Maximum GC fraction in local windows (e.g. 0.70).
    local_gc_window : int, optional
        Window size in bp for local GC (default 60).
    forbidden_sites : list of str, optional
        List of restriction sites to forbid.
    max_homopolymer : int, optional
        Maximum allowed identical consecutive nucleotides (e.g. 6).

    Returns
    -------
    str
        Optimized synonymous DNA sequence.
    """
    sequence = clean_dna(sequence)
    original_codons = get_codons(sequence)
    table = get_codon_table(organism)

    problem = pulp.LpProblem("Codon_Optimization", pulp.LpMaximize)
    variables = {}

    # 1. Decision variables: x_{pos}_{codon} in {0, 1}
    for position, original_codon in enumerate(original_codons):
        amino_acid = CODON_TO_AA[original_codon]
        synonymous_codons = table[amino_acid]
        variables[position] = {}
        for codon in synonymous_codons:
            variables[position][codon] = pulp.LpVariable(
                f"x_{position}_{codon}",
                cat=pulp.LpBinary
            )

    # 2. Constraint: Exactly one codon chosen per position
    for position in variables:
        problem += pulp.lpSum(variables[position].values()) == 1

    # 3. Objective: Maximize codon usage frequency preference
    objective_terms = []
    for position in variables:
        for codon, variable in variables[position].items():
            freq = table[CODON_TO_AA[codon]][codon]
            objective_terms.append(freq * variable)
    problem += pulp.lpSum(objective_terms)

    # 4. Global GC constraint
    total_bases = len(sequence)
    gc_expression = []
    for position in variables:
        for codon, variable in variables[position].items():
            gc_cnt = codon.count("G") + codon.count("C")
            gc_expression.append(gc_cnt * variable)
    total_gc = pulp.lpSum(gc_expression)

    if gc_min is not None:
        problem += total_gc >= gc_min * total_bases
    if gc_max is not None:
        problem += total_gc <= gc_max * total_bases

    # 5. 5' GC constraint (first W bp)
    if (five_prime_gc_min is not None or five_prime_gc_max is not None) and five_prime_window > 0:
        window_codons = min(len(original_codons), five_prime_window // 3)
        if window_codons > 0:
            window_gc_expr = []
            for position in range(window_codons):
                for codon, variable in variables[position].items():
                    gc_cnt = codon.count("G") + codon.count("C")
                    window_gc_expr.append(gc_cnt * variable)
            window_actual_bases = window_codons * 3
            if five_prime_gc_min is not None:
                problem += pulp.lpSum(window_gc_expr) >= five_prime_gc_min * window_actual_bases
            if five_prime_gc_max is not None:
                problem += pulp.lpSum(window_gc_expr) <= five_prime_gc_max * window_actual_bases

    # 6. Local sliding-window GC constraint
    if (local_gc_min is not None or local_gc_max is not None) and local_gc_window > 0:
        w_codons = local_gc_window // 3
        if 0 < w_codons <= len(original_codons):
            w_bases = w_codons * 3
            for start_idx in range(len(original_codons) - w_codons + 1):
                loc_expr = []
                for p in range(start_idx, start_idx + w_codons):
                    for codon, variable in variables[p].items():
                        gc_cnt = codon.count("G") + codon.count("C")
                        loc_expr.append(gc_cnt * variable)
                if local_gc_min is not None:
                    problem += pulp.lpSum(loc_expr) >= local_gc_min * w_bases
                if local_gc_max is not None:
                    problem += pulp.lpSum(loc_expr) <= local_gc_max * w_bases

    # 7. Forbidden restriction sites cuts
    if forbidden_sites:
        for site in forbidden_sites:
            site = site.upper().strip()
            if not site or len(site) > total_bases:
                continue
            site_len = len(site)
            for start_nt in range(total_bases - site_len + 1):
                c_start = start_nt // 3
                c_end = (start_nt + site_len - 1) // 3
                span = list(range(c_start, c_end + 1))
                if len(span) > 3:
                    continue  # Keep problem tractable for combinatorial cuts
                offset = start_nt % 3
                candidate_lists = [list(variables[p].keys()) for p in span]
                for combo in itertools.product(*candidate_lists):
                    joined = "".join(combo)
                    if joined[offset : offset + site_len] == site:
                        problem += (
                            pulp.lpSum(variables[p][combo[idx]] for idx, p in enumerate(span))
                            <= len(span) - 1
                        )

    # 8. Homopolymer cuts
    if max_homopolymer is not None and max_homopolymer > 0:
        bad_len = max_homopolymer + 1
        for base in "ACGT":
            bad_run = base * bad_len
            for start_nt in range(total_bases - bad_len + 1):
                c_start = start_nt // 3
                c_end = (start_nt + bad_len - 1) // 3
                span = list(range(c_start, c_end + 1))
                if len(span) > 3:
                    continue
                offset = start_nt % 3
                candidate_lists = [list(variables[p].keys()) for p in span]
                for combo in itertools.product(*candidate_lists):
                    joined = "".join(combo)
                    if joined[offset : offset + bad_len] == bad_run:
                        problem += (
                            pulp.lpSum(variables[p][combo[idx]] for idx, p in enumerate(span))
                            <= len(span) - 1
                        )

    # 9. Solve with CBC
    solver = pulp.PULP_CBC_CMD(msg=False)
    status = problem.solve(solver)

    if pulp.LpStatus[status] != "Optimal":
        raise RuntimeError(
            f"MILP optimization did not find an optimal solution. Solver status: {pulp.LpStatus[status]}"
        )

    # 10. Extract solution
    optimized_codons = []
    for position in range(len(original_codons)):
        selected_codon = None
        for codon, variable in variables[position].items():
            if pulp.value(variable) is not None and pulp.value(variable) > 0.5:
                selected_codon = codon
                break
        if selected_codon is None:
            raise RuntimeError(f"No codon selected at position {position}.")
        optimized_codons.append(selected_codon)

    return "".join(optimized_codons)


def milp_optimization(sequence, organism, **kwargs):
    """Alias for optimize_milp()."""
    return optimize_milp(sequence=sequence, organism=organism, **kwargs)