"""
MILP-based codon optimization.

The optimizer selects one synonymous codon for every
amino-acid position while maximizing codon usage preference.

Supported linear constraints:
- Global GC-content bounds
- 5' region GC-content bounds (exact nucleotide window)
- Local sliding-window GC-content bounds (exact nucleotide window)
- Forbidden restriction site elimination cuts (arbitrary spans and reverse complements)
- Excessive homopolymer run elimination cuts (arbitrary spans)
- Repeated subsequence elimination via lazy cut generation
"""

from typing import List, Optional, Union, Dict, Any
import pulp

from core.sequence import get_codons, clean_dna, reverse_complement
from core.verification import CODON_TO_AA
from constraints.repeats import find_repeats
from data.codon_tables import get_codon_table


def optimize_milp(
    sequence: str,
    organism: str,
    gc_min: Optional[float] = None,
    gc_max: Optional[float] = None,
    five_prime_gc_min: Optional[float] = None,
    five_prime_gc_max: Optional[float] = None,
    five_prime_window: int = 60,
    local_gc_min: Optional[float] = None,
    local_gc_max: Optional[float] = None,
    local_gc_window: int = 60,
    forbidden_sites: Optional[Union[List[str], Dict[str, str]]] = None,
    max_homopolymer: Optional[int] = None,
    repeat_length: Optional[int] = None,
    check_reverse_complement: bool = True,
    max_repeat_iterations: int = 10,
) -> str:
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
    forbidden_sites : list of str or dict, optional
        List or dict of restriction sites to forbid.
    max_homopolymer : int, optional
        Maximum allowed identical consecutive nucleotides (e.g. 6).
    repeat_length : int, optional
        Minimum length of repeated substrings to forbid via lazy cuts (e.g. 12).
    check_reverse_complement : bool, default True
        If True, also forbid reverse complements of non-palindromic restriction sites.
    max_repeat_iterations : int, default 10
        Maximum iterations of lazy repeat cut generation.

    Returns
    -------
    str
        Optimized synonymous DNA sequence.
    """
    sequence = clean_dna(sequence)
    original_codons = get_codons(sequence)
    table = get_codon_table(organism)
    total_bases = len(sequence)

    problem = pulp.LpProblem("Codon_Optimization", pulp.LpMaximize)
    variables: Dict[int, Dict[str, pulp.LpVariable]] = {}

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
    gc_expression = []
    for position in variables:
        for codon, variable in variables[position].items():
            gc_cnt = codon.count("G") + codon.count("C")
            if gc_cnt > 0:
                gc_expression.append(gc_cnt * variable)
    total_gc = pulp.lpSum(gc_expression)

    if gc_min is not None:
        problem += total_gc >= gc_min * total_bases
    if gc_max is not None:
        problem += total_gc <= gc_max * total_bases

    # 5. 5' GC constraint (first W bp, exact nucleotide window)
    if (five_prime_gc_min is not None or five_prime_gc_max is not None) and five_prime_window > 0:
        w5_bases = min(total_bases, five_prime_window)
        w5_codons = (w5_bases - 1) // 3 + 1
        window_gc_expr = []
        for position in range(w5_codons):
            j_end = min(w5_bases, 3 * position + 3)
            k1 = j_end - 3 * position
            for codon, variable in variables[position].items():
                gc_cnt = codon[:k1].count("G") + codon[:k1].count("C")
                if gc_cnt > 0:
                    window_gc_expr.append(gc_cnt * variable)
        total_win5_gc = pulp.lpSum(window_gc_expr)
        if five_prime_gc_min is not None:
            problem += total_win5_gc >= five_prime_gc_min * w5_bases
        if five_prime_gc_max is not None:
            problem += total_win5_gc <= five_prime_gc_max * w5_bases

    # 6. Local sliding-window GC constraint (exact nucleotide-level sliding window)
    if (local_gc_min is not None or local_gc_max is not None) and local_gc_window > 0:
        if total_bases >= local_gc_window:
            w_bases = local_gc_window
            for start_nt in range(total_bases - w_bases + 1):
                c_start = start_nt // 3
                c_end = (start_nt + w_bases - 1) // 3
                loc_expr = []
                for p in range(c_start, c_end + 1):
                    j_start = max(start_nt, 3 * p)
                    j_end = min(start_nt + w_bases, 3 * p + 3)
                    k0 = j_start - 3 * p
                    k1 = j_end - 3 * p
                    for codon, variable in variables[p].items():
                        gc_cnt = codon[k0:k1].count("G") + codon[k0:k1].count("C")
                        if gc_cnt > 0:
                            loc_expr.append(gc_cnt * variable)
                win_gc = pulp.lpSum(loc_expr)
                if local_gc_min is not None:
                    problem += win_gc >= local_gc_min * w_bases
                if local_gc_max is not None:
                    problem += win_gc <= local_gc_max * w_bases
        else:
            if local_gc_min is not None:
                problem += total_gc >= local_gc_min * total_bases
            if local_gc_max is not None:
                problem += total_gc <= local_gc_max * total_bases

    # Helper: Add linear cut forbidding a specific nucleotide pattern across arbitrary codon spans
    def add_forbidden_pattern_cuts(pattern: str):
        pat = pattern.upper().strip()
        pat_len = len(pat)
        if not pat or pat_len > total_bases:
            return
        for start_nt in range(total_bases - pat_len + 1):
            c_start = start_nt // 3
            c_end = (start_nt + pat_len - 1) // 3
            span = list(range(c_start, c_end + 1))
            possible = True
            cut_vars = []
            for p in span:
                j_start = max(start_nt, 3 * p)
                j_end = min(start_nt + pat_len, 3 * p + 3)
                k0 = j_start - 3 * p
                k1 = j_end - 3 * p
                sub_pat = pat[j_start - start_nt : j_end - start_nt]
                matching = [v for c, v in variables[p].items() if c[k0:k1] == sub_pat]
                if not matching:
                    possible = False
                    break
                cut_vars.extend(matching)
            if possible:
                problem.addConstraint(pulp.lpSum(cut_vars) <= len(span) - 1)

    # 7. Forbidden restriction sites cuts
    if forbidden_sites:
        sites_to_cut = set()
        if isinstance(forbidden_sites, dict):
            site_list = list(forbidden_sites.values())
        else:
            site_list = list(forbidden_sites)
        for s in site_list:
            s_clean = s.upper().strip()
            if s_clean:
                sites_to_cut.add(s_clean)
                if check_reverse_complement:
                    try:
                        sites_to_cut.add(reverse_complement(s_clean))
                    except Exception:
                        pass
        for site in sites_to_cut:
            add_forbidden_pattern_cuts(site)

    # 8. Homopolymer cuts
    if max_homopolymer is not None and max_homopolymer > 0:
        bad_len = max_homopolymer + 1
        for base in "ACGT":
            add_forbidden_pattern_cuts(base * bad_len)

    def extract_solution() -> str:
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

    # 9. Solve with CBC (with lazy repeat cuts if repeat_length is specified)
    solver = pulp.PULP_CBC_CMD(msg=False)
    max_iters = max_repeat_iterations if (repeat_length is not None and repeat_length > 0) else 1

    candidate_sequence = ""
    for iteration in range(max_iters):
        status = problem.solve(solver)
        if pulp.LpStatus[status] != "Optimal":
            raise RuntimeError(
                f"MILP optimization did not find an optimal solution. Solver status: {pulp.LpStatus[status]}"
            )
        candidate_sequence = extract_solution()
        if repeat_length is None or repeat_length <= 0:
            break

        repeats = find_repeats(candidate_sequence, repeat_length)
        if not repeats:
            break

        for rep in repeats:
            frag = rep["sequence"]
            L = len(frag)
            s1 = rep["first_position"]
            s2 = rep["second_position"]
            span1 = list(range(s1 // 3, (s1 + L - 1) // 3 + 1))
            span2 = list(range(s2 // 3, (s2 + L - 1) // 3 + 1))

            cut1 = []
            for p in span1:
                j_start = max(s1, 3 * p)
                j_end = min(s1 + L, 3 * p + 3)
                sub = frag[j_start - s1 : j_end - s1]
                cut1.extend([v for c, v in variables[p].items() if c[j_start - 3 * p : j_end - 3 * p] == sub])

            cut2 = []
            for p in span2:
                j_start = max(s2, 3 * p)
                j_end = min(s2 + L, 3 * p + 3)
                sub = frag[j_start - s2 : j_end - s2]
                cut2.extend([v for c, v in variables[p].items() if c[j_start - 3 * p : j_end - 3 * p] == sub])

            problem += pulp.lpSum(cut1) + pulp.lpSum(cut2) <= len(span1) + len(span2) - 1

    return candidate_sequence


def milp_optimization(sequence: str, organism: str, **kwargs) -> str:
    """Alias for optimize_milp()."""
    return optimize_milp(sequence=sequence, organism=organism, **kwargs)